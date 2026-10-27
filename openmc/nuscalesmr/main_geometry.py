import os
import openmc

# =====================================================================
# Sabitler (cm)
# =====================================================================
ACTIVE_HEIGHT = 200.0        # NuScale aktif yakıt boyu ~2.0 m
AXIAL_REFLECTOR = 20.0       # alt/üst su yansıtıcı kalınlığı
RADIAL_VESSEL_R = 110.0      # 7x7 kutunun köşegen yarıçapı 106.4 cm -> 110 yeterli

PIN_PITCH = 1.26
N_PINS = 17                  # 17x17 demet
ASSEMBLY_PITCH = 21.5        # 17*1.26 = 21.42 + 0.08 cm su aralığı
N_CORE = 7                   # 7x7 kafes, içinde 37 aktif demet

FUEL_R, GAP_R, CLAD_R = 0.41, 0.42, 0.47
GT_IN_R, GT_OUT_R = 0.57, 0.61          # kılavuz tüp iç/dış yarıçap
ROD_ABS_R, ROD_CLAD_R = 0.38, 0.42      # kontrol çubuğu emici / kaplama

# Standart 17x17 düzeninde 24 kılavuz tüp + merkezde 1 enstrüman tüpü
GUIDE_TUBES = [
    (2, 5), (2, 8), (2, 11), (3, 3), (3, 13),
    (5, 2), (5, 5), (5, 8), (5, 11), (5, 14),
    (8, 2), (8, 5), (8, 11), (8, 14),
    (11, 2), (11, 5), (11, 8), (11, 11), (11, 14),
    (13, 3), (13, 13), (14, 5), (14, 8), (14, 11),
]
INSTRUMENT_TUBE = (8, 8)

# 37 aktif demet: satır -> (ilk_kolon, son_kolon)
ROW_SPAN = {0: (2, 4), 1: (1, 5), 2: (0, 6), 3: (0, 6), 4: (0, 6), 5: (1, 5), 6: (2, 4)}

# 16 kontrol çubuklu demet (4 bank x 4 demet, 90 derece dönme simetrisi).
# NOT: Gerçek NuScale CRA haritanız varsa burayı onunla değiştirin.
BANK_MAP = {
    (2, 1): 0, (5, 2): 0, (4, 5): 0, (1, 4): 0,   # bank 1
    (2, 5): 1, (1, 2): 1, (4, 1): 1, (5, 4): 1,   # bank 2
    (3, 5): 2, (1, 3): 2, (3, 1): 2, (5, 3): 2,   # bank 3
    (5, 5): 3, (1, 5): 3, (1, 1): 3, (5, 1): 3,   # bank 4
}


# =====================================================================
# Malzemeler
# =====================================================================
def build_materials(coolant_density=0.70):
    fuel = openmc.Material(name='fuel')
    fuel.add_nuclide('U235', 0.045)
    fuel.add_nuclide('U238', 0.955)
    fuel.add_nuclide('O16', 2.0)
    fuel.set_density('g/cm3', 10.5)

    gap = openmc.Material(name='gap')          # su değil, helyum dolgu
    gap.add_nuclide('He4', 1.0)
    gap.set_density('g/cm3', 0.0016)

    clad = openmc.Material(name='cladding')
    clad.add_element('Zr', 1.0)
    clad.set_density('g/cm3', 6.5)

    steel = openmc.Material(name='ss304')
    steel.add_element('Fe', 0.70, percent_type='wo')
    steel.add_element('Cr', 0.19, percent_type='wo')
    steel.add_element('Ni', 0.11, percent_type='wo')
    steel.set_density('g/cm3', 8.0)

    # B4C: 4 B : 1 C  (orijinal kodda oranlar ters yazılmıştı)
    control = openmc.Material(name='control_B4C')
    control.add_element('B', 0.8)
    control.add_nuclide('C12', 0.2)
    control.set_density('g/cm3', 1.7)

    water = openmc.Material(name='moderator')
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.set_density('g/cm3', coolant_density)
    water.add_s_alpha_beta('c_H_in_H2O')

    refl = openmc.Material(name='reflector')
    refl.add_nuclide('H1', 2.0)
    refl.add_nuclide('O16', 1.0)
    refl.set_density('g/cm3', 1.0)
    refl.add_s_alpha_beta('c_H_in_H2O')

    return dict(fuel=fuel, gap=gap, clad=clad, steel=steel,
                control=control, water=water, refl=refl)


# =====================================================================
# Geometri
# =====================================================================
def build_model(bank_insertion=(0.0, 0.0, 0.0, 0.0), coolant_density=0.70):
    """bank_insertion: 4 bank için 0 (tam çekili) .. 1 (tam inik) eklem oranı."""
    m = build_materials(coolant_density)
    z_a0, z_a1 = -ACTIVE_HEIGHT / 2.0, ACTIVE_HEIGHT / 2.0

    # ---- yakıt çubuğu ------------------------------------------------
    s_fuel = openmc.ZCylinder(r=FUEL_R)
    s_gap = openmc.ZCylinder(r=GAP_R)
    s_clad = openmc.ZCylinder(r=CLAD_R)
    pin = openmc.Universe(name='fuel_pin', cells=[
        openmc.Cell(fill=m['fuel'], region=-s_fuel),
        openmc.Cell(fill=m['gap'], region=+s_fuel & -s_gap),
        openmc.Cell(fill=m['clad'], region=+s_gap & -s_clad),
        openmc.Cell(fill=m['water'], region=+s_clad),
    ])

    # ---- kılavuz tüp (boş) ve kontrol çubuklu kılavuz tüp -----------------
    s_gt_in = openmc.ZCylinder(r=GT_IN_R)
    s_gt_out = openmc.ZCylinder(r=GT_OUT_R)
    s_abs = openmc.ZCylinder(r=ROD_ABS_R)
    s_rclad = openmc.ZCylinder(r=ROD_CLAD_R)

    gt_empty = openmc.Universe(name='guide_tube_empty', cells=[
        openmc.Cell(fill=m['water'], region=-s_gt_in),
        openmc.Cell(fill=m['clad'], region=+s_gt_in & -s_gt_out),
        openmc.Cell(fill=m['water'], region=+s_gt_out),
    ])
    gt_rod = openmc.Universe(name='guide_tube_rod', cells=[
        openmc.Cell(fill=m['control'], region=-s_abs),
        openmc.Cell(fill=m['steel'], region=+s_abs & -s_rclad),
        openmc.Cell(fill=m['water'], region=+s_rclad & -s_gt_in),
        openmc.Cell(fill=m['clad'], region=+s_gt_in & -s_gt_out),
        openmc.Cell(fill=m['water'], region=+s_gt_out),
    ])

    def rod_axial(insertion):
        """Çubuk üstten girer; ucu z_a1 - insertion*H konumunda."""
        insertion = min(max(float(insertion), 0.0), 1.0)
        if insertion <= 0.0:
            return gt_empty
        if insertion >= 1.0:
            return gt_rod
        zp = openmc.ZPlane(z0=z_a1 - insertion * ACTIVE_HEIGHT)
        return openmc.Universe(cells=[
            openmc.Cell(fill=gt_rod, region=+zp),
            openmc.Cell(fill=gt_empty, region=-zp),
        ])

    # ---- sınırsız su evreni (kafes 'outer' için) -------------------------
    water_univ = openmc.Universe(cells=[openmc.Cell(fill=m['water'])])
    refl_univ = openmc.Universe(cells=[openmc.Cell(fill=m['refl'])])

    # ---- 17x17 demet -------------------------------------------------------
    def make_assembly(gt_univ):
        layout = [[pin for _ in range(N_PINS)] for _ in range(N_PINS)]
        for i, j in GUIDE_TUBES:
            layout[i][j] = gt_univ
        i, j = INSTRUMENT_TUBE
        layout[i][j] = gt_empty

        lat = openmc.RectLattice()
        lat.pitch = (PIN_PITCH, PIN_PITCH)
        lat.lower_left = (-N_PINS * PIN_PITCH / 2.0, -N_PINS * PIN_PITCH / 2.0)
        lat.universes = layout
        lat.outer = water_univ      # 21.42 -> 21.5 cm arasındaki su aralığı
        # Hücrede bölge yok: evren tüm uzayı kaplar, sınırı üst kafes belirler.
        return openmc.Universe(cells=[openmc.Cell(fill=lat)])

    asm_plain = make_assembly(gt_empty)
    asm_bank = [make_assembly(rod_axial(x)) for x in bank_insertion]

    # ---- çekirdek kafesi (7x7, 37 aktif demet) -----------------------------
    core_u = [[refl_univ for _ in range(N_CORE)] for _ in range(N_CORE)]
    n_active = 0
    for i, (a, b) in ROW_SPAN.items():
        for j in range(a, b + 1):
            core_u[i][j] = asm_bank[BANK_MAP[(i, j)]] if (i, j) in BANK_MAP else asm_plain
            n_active += 1
    assert n_active == 37, n_active
    assert len(BANK_MAP) == 16

    half = N_CORE * ASSEMBLY_PITCH / 2.0
    core_lat = openmc.RectLattice()
    core_lat.pitch = (ASSEMBLY_PITCH, ASSEMBLY_PITCH)
    core_lat.lower_left = (-half, -half)
    core_lat.universes = core_u
    core_lat.outer = refl_univ

    # ---- dış sınırlar ------------------------------------------------------
    bx0, bx1 = openmc.XPlane(x0=-half), openmc.XPlane(x0=half)
    by0, by1 = openmc.YPlane(y0=-half), openmc.YPlane(y0=half)
    box = +bx0 & -bx1 & +by0 & -by1

    vessel = openmc.ZCylinder(r=RADIAL_VESSEL_R, boundary_type='vacuum')
    z_bot = openmc.ZPlane(z0=z_a0 - AXIAL_REFLECTOR, boundary_type='vacuum')
    z_top = openmc.ZPlane(z0=z_a1 + AXIAL_REFLECTOR, boundary_type='vacuum')
    z0p = openmc.ZPlane(z0=z_a0)
    z1p = openmc.ZPlane(z0=z_a1)

    core_cell = openmc.Cell(name='core', fill=core_lat, region=box & +z0p & -z1p)
    radial_refl = openmc.Cell(name='radial_reflector', fill=m['refl'],
                              region=-vessel & ~box & +z0p & -z1p)
    bottom_refl = openmc.Cell(name='bottom_reflector', fill=m['refl'],
                              region=-vessel & +z_bot & -z0p)
    top_refl = openmc.Cell(name='top_reflector', fill=m['refl'],
                           region=-vessel & +z1p & -z_top)

    root = openmc.Universe(cells=[core_cell, radial_refl, bottom_refl, top_refl])
    geometry = openmc.Geometry(root)
    materials = openmc.Materials(list(m.values()))
    return geometry, materials


def export(geometry, materials, outdir='geometry'):
    os.makedirs(outdir, exist_ok=True)
    geometry.export_to_xml(os.path.join(outdir, 'geometry.xml'))
    # materials.xml AYNI Material nesnelerinden üretilmeli (ID uyumu için)
    materials.export_to_xml(os.path.join(outdir, 'materials.xml'))


if __name__ == '__main__':
    geom, mats = build_model(bank_insertion=(0.0, 0.0, 0.0, 0.0), coolant_density=0.70)
    export(geom, mats)