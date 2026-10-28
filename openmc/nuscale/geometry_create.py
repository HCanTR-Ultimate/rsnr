import openmc
import os

ACTIVE_HEIGHT = 200.0
N_PINS = 17
PIN_PITCH = 1.26
ASSEMBLY_PITCH = 21.5

R_FUEL = 0.4096
R_CLAD_IN = 0.4180
R_CLAD_OUT = 0.4750
R_GT_IN = 0.5715
R_GT_OUT = 0.6121
R_ABS = 0.4331
R_ABS_CLAD_IN = 0.4369
R_ABS_CLAD_OUT = 0.4839

Z_BOTTOM = -60.0
Z_LOWER_PLATE = (-30.0, -10.0)
Z_ENDPLUG = (-10.0, 0.0)
Z_PLENUM = (ACTIVE_HEIGHT, ACTIVE_HEIGHT + 15.0)
Z_UPPER_PLATE = (215.0, 235.0)
Z_TOP = 265.0

R_REFLECTOR = 90.0
R_BARREL_OUT = 92.5
R_VESSEL_IN = 132.0
R_VESSEL_OUT = 137.2

T_FUEL = 900.0
T_CLAD = 600.0
T_STRUCT = 550.0
T_COOLANT = 544.15

ENRICHMENT = 4.55
BORON_PPM = 1000.0
COOLANT_DENSITY = 0.772
RHO_UO2 = 10.40
RHO_UO2_GD = 10.20
GD_ENRICHMENT = 2.0
GD_WT = 0.08

SS304 = ({'Fe': 0.695, 'Cr': 0.19, 'Ni': 0.095, 'Mn': 0.02}, 7.9)
CARBON_STEEL = ({'Fe': 0.9675, 'Mn': 0.014, 'Ni': 0.007, 'Mo': 0.005, 'Cr': 0.003, 'C': 0.0035}, 7.85)
M5 = ({'Zr': 0.99, 'Nb': 0.01}, 6.55)
AIC = ({'Ag': 0.80, 'In': 0.15, 'Cd': 0.05}, 10.17)
HELIUM = ({'He': 1.0}, 0.0020)

GUIDE_TUBES = {
    (2, 5), (2, 8), (2, 11), (3, 3), (3, 13), (5, 2), (5, 5), (5, 8), (5, 11),
    (5, 14), (8, 2), (8, 5), (8, 11), (8, 14), (11, 2), (11, 5), (11, 8),
    (11, 11), (11, 14), (13, 3), (13, 13), (14, 5), (14, 8), (14, 11),
}
INSTRUMENT_TUBE = (8, 8)
GD_PINS = {
    (1, 5), (1, 11), (5, 1), (11, 1), (15, 5), (15, 11), (5, 15), (11, 15),
    (2, 2), (2, 14), (14, 2), (14, 14), (4, 4), (4, 12), (12, 4), (12, 12),
}

ASSEMBLIES = [
    (dr, dc)
    for dr in range(-3, 4)
    for dc in range(-3, 4)
    if (abs(dr) <= 1 and abs(dc) <= 3)
    or (abs(dr) == 2 and abs(dc) <= 2)
    or (abs(dr) == 3 and abs(dc) <= 1)
]


def _orbit(a, b):
    return [(a, b), (-b, a), (-a, -b), (b, -a)]


CRA_POSITIONS = {p for seed in [(0, 2), (1, 1), (1, 3), (2, 1)] for p in _orbit(*seed)}


def water_recipe(rho, ppm):
    b = ppm * 1e-6
    comp = {'H': 0.111894 * (1.0 - b), 'O': 0.888106 * (1.0 - b)}
    if b > 0.0:
        comp['B'] = b
    return comp, rho


def mix(name, parts, temperature):
    comp = {}
    total = 0.0
    for (recipe, rho), vf in parts:
        total += rho * vf
        for el, w in recipe.items():
            comp[el] = comp.get(el, 0.0) + rho * vf * w
    mat = openmc.Material(name=name)
    for el, m in comp.items():
        mat.add_element(el, m / total, 'wo')
    mat.set_density('g/cm3', total)
    if 'H' in comp:
        mat.add_s_alpha_beta('c_H_in_H2O')
    mat.temperature = temperature
    return mat


def make_uo2(name, enrichment, rho):
    mat = openmc.Material(name=name)
    mat.add_elements_from_formula('UO2', enrichment=enrichment)
    mat.set_density('g/cm3', rho)
    mat.temperature = T_FUEL
    return mat


def make_gd_fuel():
    uo2 = make_uo2('uo2_gd_base', GD_ENRICHMENT, 10.96)
    gd2o3 = openmc.Material(name='gd2o3')
    gd2o3.add_elements_from_formula('Gd2O3')
    gd2o3.set_density('g/cm3', 7.41)
    mat = openmc.Material.mix_materials(
        [uo2, gd2o3], [1.0 - GD_WT, GD_WT], 'wo', name='fuel_gd')
    mat.set_density('g/cm3', RHO_UO2_GD)
    mat.temperature = T_FUEL
    return mat


def rings(name, mats, radii, outer):
    surfs = [openmc.ZCylinder(r=r) for r in radii]
    cells = []
    for i, m in enumerate(mats):
        region = -surfs[0] if i == 0 else (+surfs[i - 1] & -surfs[i])
        cells.append(openmc.Cell(name=f'{name}_{i}', fill=m, region=region))
    cells.append(openmc.Cell(name=f'{name}_out', fill=outer, region=+surfs[-1]))
    return openmc.Universe(name=name, cells=cells)


class _Core:
    def __init__(self, enrichment, boron_ppm, coolant_density):
        self.enrichment = enrichment
        self.boron = boron_ppm
        self.rho = coolant_density
        self._planes = {}
        self._lattices = {}
        self._assemblies = {}
        self._build_materials()
        self._build_pins()
        self.reflector = self._build_reflector()
        self.geometry = self._build_geometry()

    def _plane(self, z):
        if z not in self._planes:
            self._planes[z] = openmc.ZPlane(z0=z)
        return self._planes[z]

    def _region(self, lo, hi):
        region = None
        if lo is not None:
            region = +self._plane(lo)
        if hi is not None:
            upper = -self._plane(hi)
            region = upper if region is None else region & upper
        return region

    def _build_materials(self):
        water = water_recipe(self.rho, self.boron)
        self.fuel = make_uo2('fuel', self.enrichment, RHO_UO2)
        self.fuel_gd = make_gd_fuel()
        self.coolant = mix('coolant', [(water, 1.0)], T_COOLANT)
        self.clad = mix('clad_m5', [(M5, 1.0)], T_CLAD)
        self.helium = mix('helium', [(HELIUM, 1.0)], T_CLAD)
        self.steel = mix('ss304', [(SS304, 1.0)], T_STRUCT)
        self.vessel_steel = mix('vessel_steel', [(CARBON_STEEL, 1.0)], T_STRUCT)
        self.aic = mix('aic', [(AIC, 1.0)], T_CLAD)
        self.plate = mix('core_plate', [(SS304, 0.30), (water, 0.70)], T_STRUCT)
        self.reflector_mat = mix('reflector', [(SS304, 0.80), (water, 0.20)], T_STRUCT)

    def _build_pins(self):
        cool = self.coolant
        self.water_u = openmc.Universe(cells=[openmc.Cell(fill=cool)])
        self.fuel_u = rings(
            'fuel_pin', [self.fuel, self.helium, self.clad],
            [R_FUEL, R_CLAD_IN, R_CLAD_OUT], cool)
        self.gd_u = rings(
            'gd_pin', [self.fuel_gd, self.helium, self.clad],
            [R_FUEL, R_CLAD_IN, R_CLAD_OUT], cool)
        self.gt_u = rings(
            'guide_tube', [cool, self.clad], [R_GT_IN, R_GT_OUT], cool)
        self.rod_u = rings(
            'control_rod',
            [self.aic, self.helium, self.steel, cool, self.clad],
            [R_ABS, R_ABS_CLAD_IN, R_ABS_CLAD_OUT, R_GT_IN, R_GT_OUT], cool)
        self.endplug_u = rings('end_plug', [self.clad], [R_CLAD_OUT], cool)
        self.plenum_u = rings(
            'plenum', [self.helium, self.clad], [R_CLAD_IN, R_CLAD_OUT], cool)

    def _lattice(self, kind, rod_in):
        key = (kind, rod_in)
        if key in self._lattices:
            return self._lattices[key]
        rows = []
        for i in range(N_PINS):
            row = []
            for j in range(N_PINS):
                p = (i, j)
                if p in GUIDE_TUBES:
                    u = self.rod_u if rod_in else self.gt_u
                elif p == INSTRUMENT_TUBE:
                    u = self.gt_u
                elif kind == 'active':
                    u = self.gd_u if p in GD_PINS else self.fuel_u
                elif kind == 'plenum':
                    u = self.plenum_u
                else:
                    u = self.endplug_u
                row.append(u)
            rows.append(row)
        lat = openmc.RectLattice(name=f'{kind}_{"rodded" if rod_in else "plain"}')
        lat.lower_left = (-0.5 * N_PINS * PIN_PITCH, -0.5 * N_PINS * PIN_PITCH)
        lat.pitch = (PIN_PITCH, PIN_PITCH)
        lat.universes = rows
        lat.outer = self.water_u
        self._lattices[key] = lat
        return lat

    def _assembly(self, has_cra):
        if has_cra in self._assemblies:
            return self._assemblies[has_cra]
        segs = [
            (None, Z_LOWER_PLATE[0], self.coolant),
            (Z_LOWER_PLATE[0], Z_LOWER_PLATE[1], self.plate),
            (Z_ENDPLUG[0], Z_ENDPLUG[1], self._lattice('endplug', False)),
            (0.0, ACTIVE_HEIGHT, self._lattice('active', False)),
            (Z_PLENUM[0], Z_PLENUM[1], self._lattice('plenum', has_cra)),
            (Z_UPPER_PLATE[0], Z_UPPER_PLATE[1], self.plate),
            (Z_UPPER_PLATE[1], None, self.coolant),
        ]
        cells = [openmc.Cell(fill=f, region=self._region(lo, hi)) for lo, hi, f in segs]
        uni = openmc.Universe(name='assembly_cra' if has_cra else 'assembly', cells=cells)
        self._assemblies[has_cra] = uni
        return uni

    def _build_reflector(self):
        segs = [
            (None, Z_LOWER_PLATE[0], self.coolant),
            (Z_LOWER_PLATE[0], Z_LOWER_PLATE[1], self.plate),
            (Z_ENDPLUG[0], Z_UPPER_PLATE[0], self.reflector_mat),
            (Z_UPPER_PLATE[0], Z_UPPER_PLATE[1], self.plate),
            (Z_UPPER_PLATE[1], None, self.coolant),
        ]
        cells = [openmc.Cell(fill=m, region=self._region(lo, hi)) for lo, hi, m in segs]
        return openmc.Universe(name='reflector', cells=cells)

    def _build_geometry(self):
        assembly_set = set(ASSEMBLIES)
        rows = []
        for dr in range(-3, 4):
            row = []
            for dc in range(-3, 4):
                if (dr, dc) in assembly_set:
                    row.append(self._assembly((dr, dc) in CRA_POSITIONS))
                else:
                    row.append(self.reflector)
            rows.append(row)
        core = openmc.RectLattice(name='core')
        core.lower_left = (-3.5 * ASSEMBLY_PITCH, -3.5 * ASSEMBLY_PITCH)
        core.pitch = (ASSEMBLY_PITCH, ASSEMBLY_PITCH)
        core.universes = rows
        core.outer = self.reflector

        zmin = openmc.ZPlane(z0=Z_BOTTOM, boundary_type='vacuum')
        zmax = openmc.ZPlane(z0=Z_TOP, boundary_type='vacuum')
        r_ref = openmc.ZCylinder(r=R_REFLECTOR)
        r_barrel = openmc.ZCylinder(r=R_BARREL_OUT)
        r_vin = openmc.ZCylinder(r=R_VESSEL_IN)
        r_vout = openmc.ZCylinder(r=R_VESSEL_OUT, boundary_type='vacuum')
        axial = +zmin & -zmax
        cells = [
            openmc.Cell(name='core', fill=core, region=-r_ref & axial),
            openmc.Cell(name='barrel', fill=self.steel, region=+r_ref & -r_barrel & axial),
            openmc.Cell(name='downcomer', fill=self.coolant, region=+r_barrel & -r_vin & axial),
            openmc.Cell(name='vessel', fill=self.vessel_steel, region=+r_vin & -r_vout & axial),
        ]
        return openmc.Geometry(cells)


def build_geometry(enrichment=ENRICHMENT, boron_ppm=BORON_PPM,
                   coolant_density=COOLANT_DENSITY):
    core = _Core(enrichment, boron_ppm, coolant_density)
    materials = openmc.Materials(list(core.geometry.get_all_materials().values()))
    return core.geometry, materials


if __name__ == '__main__':
    os.makedirs('resources', exist_ok=True)
    geometry, materials = build_geometry()
    geometry.export_to_xml('resources/geometry.xml')
    materials.export_to_xml('resources/materials.xml')