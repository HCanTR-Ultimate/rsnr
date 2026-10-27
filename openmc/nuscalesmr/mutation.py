import json
import math
import openmc
import os
import openmc.deplete as dep
from physics import (get_water_density, get_water_temperature_k, get_axial_density_profile,
                     make_fuel_material, fresh_fuel_material, atom_densities, EXCLUDED, TABLE_PATH)

N_AXIAL = 6
ACTIVE_HEIGHT = 200.0
AXIAL_REFLECTOR = 20.0
RADIAL_VESSEL_R = 110.0
PIN_PITCH = 1.26
N_PINS = 17
ASSEMBLY_PITCH = 21.5
N_CORE = 7
FUEL_R, GAP_R, CLAD_R = 0.41, 0.42, 0.47
GT_IN_R, GT_OUT_R = 0.57, 0.61
ROD_ABS_R, ROD_CLAD_R = 0.38, 0.42

FUEL_TEMP_K = 900.0
CLAD_DT_K = 25.0
T_REFLECTOR_C = 258.0
BORON_PPM = 0.0
TEMPERATURE_SETTINGS = {'method': 'interpolation', 'default': 294.0, 'range': (250.0, 1300.0)}

CORE_POWER_MWT = 250.0
T_COOLANT_AVG_C = 285.0
DEPLETION_STEPS = [0.5, 0.5, 1.0, 2.0, 3.0, 3.0] + [3.0] * 10

GUIDE_TUBES = [
    (2, 5), (2, 8), (2, 11), (3, 3), (3, 13),
    (5, 2), (5, 5), (5, 8), (5, 11), (5, 14),
    (8, 2), (8, 5), (8, 11), (8, 14),
    (11, 2), (11, 5), (11, 8), (11, 11), (11, 14),
    (13, 3), (13, 13), (14, 5), (14, 8), (14, 11),
]
INSTRUMENT_TUBE = (8, 8)
ROW_SPAN = {0: (2, 4), 1: (1, 5), 2: (0, 6), 3: (0, 6), 4: (0, 6), 5: (1, 5), 6: (2, 4)}
BANK_MAP = {
    (2, 1): 0, (5, 2): 0, (4, 5): 0, (1, 4): 0,
    (2, 5): 1, (1, 2): 1, (4, 1): 1, (5, 4): 1,
    (3, 5): 2, (1, 3): 2, (3, 1): 2, (5, 3): 2,
    (5, 5): 3, (1, 5): 3, (1, 1): 3, (5, 1): 3,
}


def make_water(name, density, boron_ppm=BORON_PPM):
    w = openmc.Material(name=name)
    w.add_nuclide('H1', 2.0)
    w.add_nuclide('O16', 1.0)
    if boron_ppm > 0.0:
        wb = boron_ppm * 1.0e-6
        w.add_element('B', (wb / 10.811) / ((1.0 - wb) / 18.015))
    w.set_density('g/cm3', density)
    w.add_s_alpha_beta('c_H_in_H2O')
    return w


def make_shared_materials():
    gap = openmc.Material(name='gap_He')
    gap.add_nuclide('He4', 1.0)
    gap.set_density('g/cm3', 0.0016)

    clad = openmc.Material(name='cladding_M5')
    clad.add_element('Zr', 0.99, percent_type='wo')
    clad.add_element('Nb', 0.01, percent_type='wo')
    clad.set_density('g/cm3', 6.5)

    steel = openmc.Material(name='ss304')
    steel.add_element('Fe', 0.70, percent_type='wo')
    steel.add_element('Cr', 0.19, percent_type='wo')
    steel.add_element('Ni', 0.11, percent_type='wo')
    steel.set_density('g/cm3', 8.0)

    control = openmc.Material(name='control_B4C')
    control.add_element('B', 0.8)
    control.add_nuclide('C12', 0.2)
    control.set_density('g/cm3', 1.7)

    refl = make_water('reflector', get_water_density(T_REFLECTOR_C))
    return dict(gap=gap, clad=clad, steel=steel, control=control, refl=refl)


def _cell(fill, region=None, T=None, name=''):
    c = openmc.Cell(fill=fill, region=region, name=name)
    if T is not None:
        c.temperature = T
    return c


def build_model(control_rods, xe_avg, xe_gradient, coolant_density_avg,
                coolant_gradient, burnup, n_axial=N_AXIAL):
    if len(control_rods) != 4:
        raise ValueError("control_rods 4 bank için 4 değer içermeli")

    rho = get_axial_density_profile(coolant_density_avg, coolant_gradient, n_axial)
    xe = get_axial_density_profile(xe_avg, xe_gradient, n_axial)
    t_water = [get_water_temperature_k(r) for r in rho]

    shared = make_shared_materials()
    gap, clad, steel, control, refl = (shared[k] for k in ('gap', 'clad', 'steel', 'control', 'refl'))
    t_refl = T_REFLECTOR_C + 273.15

    zs = [-ACTIVE_HEIGHT / 2.0 + i * ACTIVE_HEIGHT / n_axial for i in range(n_axial + 1)]
    zplanes = [openmc.ZPlane(z0=z) for z in zs]
    z_a1 = zs[-1]

    s_fuel, s_gap, s_clad = (openmc.ZCylinder(r=r) for r in (FUEL_R, GAP_R, CLAD_R))
    s_gt_in, s_gt_out = openmc.ZCylinder(r=GT_IN_R), openmc.ZCylinder(r=GT_OUT_R)
    s_abs, s_rclad = openmc.ZCylinder(r=ROD_ABS_R), openmc.ZCylinder(r=ROD_CLAD_R)

    all_materials = [gap, clad, steel, control, refl]
    refl_univ = openmc.Universe(cells=[_cell(refl, None, t_refl)])

    def build_layer(k):
        fuel = make_fuel_material(f'fuel_{k}', burnup, xe[k])
        water = make_water(f'moderator_{k}', rho[k])
        all_materials.extend([fuel, water])
        t_w = t_water[k]
        t_clad = t_w + CLAD_DT_K
        t_gap = 0.5 * (FUEL_TEMP_K + t_clad)

        pin = openmc.Universe(cells=[
            _cell(fuel, -s_fuel, FUEL_TEMP_K),
            _cell(gap, +s_fuel & -s_gap, t_gap),
            _cell(clad, +s_gap & -s_clad, t_clad),
            _cell(water, +s_clad, t_w),
        ])
        gt_empty = openmc.Universe(cells=[
            _cell(water, -s_gt_in, t_w),
            _cell(clad, +s_gt_in & -s_gt_out, t_clad),
            _cell(water, +s_gt_out, t_w),
        ])
        gt_rod = openmc.Universe(cells=[
            _cell(control, -s_abs, t_w),
            _cell(steel, +s_abs & -s_rclad, t_w),
            _cell(water, +s_rclad & -s_gt_in, t_w),
            _cell(clad, +s_gt_in & -s_gt_out, t_clad),
            _cell(water, +s_gt_out, t_w),
        ])
        water_univ = openmc.Universe(cells=[_cell(water, None, t_w)])
        zk0, zk1 = zs[k], zs[k + 1]

        def rod_axial(ins):
            ins = min(max(float(ins), 0.0), 1.0)
            tip = z_a1 - ins * ACTIVE_HEIGHT
            if ins <= 0.0 or tip >= zk1:
                return gt_empty
            if tip <= zk0:
                return gt_rod
            zp = openmc.ZPlane(z0=tip)
            return openmc.Universe(cells=[openmc.Cell(fill=gt_rod, region=+zp),
                                          openmc.Cell(fill=gt_empty, region=-zp)])

        def make_assembly(gt_univ):
            layout = [[pin for _ in range(N_PINS)] for _ in range(N_PINS)]
            for i, j in GUIDE_TUBES:
                layout[i][j] = gt_univ
            layout[INSTRUMENT_TUBE[0]][INSTRUMENT_TUBE[1]] = gt_empty
            lat = openmc.RectLattice()
            lat.pitch = (PIN_PITCH, PIN_PITCH)
            lat.lower_left = (-N_PINS * PIN_PITCH / 2.0, -N_PINS * PIN_PITCH / 2.0)
            lat.universes = layout
            lat.outer = water_univ
            return openmc.Universe(cells=[openmc.Cell(fill=lat)])

        asm_plain = make_assembly(gt_empty)
        asm_bank = [make_assembly(rod_axial(x)) for x in control_rods]

        core_u = [[refl_univ for _ in range(N_CORE)] for _ in range(N_CORE)]
        for i, (a, b) in ROW_SPAN.items():
            for j in range(a, b + 1):
                core_u[i][j] = asm_bank[BANK_MAP[(i, j)]] if (i, j) in BANK_MAP else asm_plain

        half = N_CORE * ASSEMBLY_PITCH / 2.0
        lat = openmc.RectLattice()
        lat.pitch = (ASSEMBLY_PITCH, ASSEMBLY_PITCH)
        lat.lower_left = (-half, -half)
        lat.universes = core_u
        lat.outer = refl_univ
        return lat

    half = N_CORE * ASSEMBLY_PITCH / 2.0
    box = (+openmc.XPlane(x0=-half) & -openmc.XPlane(x0=half)
           & +openmc.YPlane(y0=-half) & -openmc.YPlane(y0=half))
    vessel = openmc.ZCylinder(r=RADIAL_VESSEL_R, boundary_type='vacuum')
    z_bot = openmc.ZPlane(z0=zs[0] - AXIAL_REFLECTOR, boundary_type='vacuum')
    z_top = openmc.ZPlane(z0=zs[-1] + AXIAL_REFLECTOR, boundary_type='vacuum')

    cells = [openmc.Cell(name=f'core_layer_{k}', fill=build_layer(k),
                         region=box & +zplanes[k] & -zplanes[k + 1]) for k in range(n_axial)]
    cells.append(_cell(refl, -vessel & ~box & +zplanes[0] & -zplanes[-1], t_refl, 'radial_reflector'))
    cells.append(_cell(refl, -vessel & +z_bot & -zplanes[0], t_refl, 'bottom_reflector'))
    cells.append(_cell(refl, -vessel & +zplanes[-1] & -z_top, t_refl, 'top_reflector'))

    return openmc.Materials(all_materials), openmc.Geometry(openmc.Universe(cells=cells))


def mutate_model(control_rods, xe_avg, xe_gradient, coolant_density_avg,
                 coolant_gradient, burnup):
    openmc.reset_auto_ids()
    return build_model(control_rods, xe_avg, xe_gradient,
                       coolant_density_avg, coolant_gradient, burnup)


def generate_burnup_table():
    openmc.reset_auto_ids()
    rho = get_water_density(T_COOLANT_AVG_C)
    t_w = get_water_temperature_k(rho)
    t_clad = t_w + CLAD_DT_K

    fuel = fresh_fuel_material('fuel')
    fuel.depletable = True
    fuel.volume = math.pi * FUEL_R ** 2 * ACTIVE_HEIGHT
    m_hm_core = fuel.fissionable_mass * 264 * 37
    spec_power = CORE_POWER_MWT * 1.0e6 / m_hm_core

    sh = make_shared_materials()
    water = make_water('moderator', rho)
    s_f, s_g, s_c = (openmc.ZCylinder(r=r) for r in (FUEL_R, GAP_R, CLAD_R))
    p = PIN_PITCH / 2.0
    h = ACTIVE_HEIGHT / 2.0
    prism = (+openmc.XPlane(x0=-p, boundary_type='reflective') & -openmc.XPlane(x0=p, boundary_type='reflective')
             & +openmc.YPlane(y0=-p, boundary_type='reflective') & -openmc.YPlane(y0=p, boundary_type='reflective')
             & +openmc.ZPlane(z0=-h, boundary_type='reflective') & -openmc.ZPlane(z0=h, boundary_type='reflective'))
    cells = [
        _cell(fuel, -s_f & prism, FUEL_TEMP_K),
        _cell(sh['gap'], +s_f & -s_g & prism, 0.5 * (FUEL_TEMP_K + t_clad)),
        _cell(sh['clad'], +s_g & -s_c & prism, t_clad),
        _cell(water, +s_c & prism, t_w),
    ]
    os.environ["OMP_NUM_THREADS"] = "11"
    settings = openmc.Settings()
    settings.batches, settings.inactive, settings.particles = 150, 30, 20000
    settings.temperature = TEMPERATURE_SETTINGS
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-p, -p, -h), (p, p, h), only_fissionable=True))
    model = openmc.Model(openmc.Geometry(openmc.Universe(cells=cells)),
                         openmc.Materials([fuel, sh['gap'], sh['clad'], water]), settings)

    op = dep.CoupledOperator(model)
    dep.CECMIntegrator(op, DEPLETION_STEPS, power_density=spec_power,
                       timestep_units='MWd/kg').integrate()

    try:
        res = dep.Results('depletion_results.h5')
    except TypeError:
        res = dep.Results.from_hdf5('depletion_results.h5')

    lib = openmc.data.DataLibrary.from_xml()
    have = {n for e in lib if e['type'] == 'neutron' for n in e['materials']}

    burnups = [0.0]
    for s in DEPLETION_STEPS:
        burnups.append(burnups[-1] + s)
    comps, dropped = [], set()
    for i in range(len(burnups)):
        kept = {}
        for n, d in atom_densities(res.export_to_materials(i)[0]).items():
            if n in EXCLUDED or d <= 0.0:
                continue
            if n not in have:
                dropped.add(n)
                continue
            kept[n] = d
        comps.append(kept)

    with open(TABLE_PATH, 'w') as f:
        json.dump({'burnup_gwd_t': burnups, 'specific_power_w_per_ghm': spec_power,
                   'compositions': comps}, f)
    print(f"Yazıldı: {TABLE_PATH}")
    print(f"Verisi olmadığı için atılan nüklitler: {sorted(dropped)}")


if __name__ == '__main__':
    generate_burnup_table()