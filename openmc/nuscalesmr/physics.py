import json
import math
import os
import openmc
from CoolProp.CoolProp import PropsSI

PRESSURE_PA = 13.8e6
ENRICHMENT_WO = 4.5
FUEL_DENSITY = 10.5
EXCLUDED = ('Xe135', 'I135')
TABLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'burnup_table.json')

LAMBDA_I135 = math.log(2.0) / (6.57 * 3600.0)
LAMBDA_XE135 = math.log(2.0) / (9.14 * 3600.0)
GAMMA_I135 = 0.0639
GAMMA_XE135 = 0.0025
ENERGY_PER_FISSION_J = 200.0e6 * 1.602176634e-19

_cache = {}


def get_water_density(temperature_c, pressure_pa=PRESSURE_PA):
    return PropsSI('D', 'T', temperature_c + 273.15, 'P', pressure_pa, 'Water') / 1000.0


def get_water_temperature_k(density_g_cm3, pressure_pa=PRESSURE_PA):
    rho_sat = PropsSI('D', 'P', pressure_pa, 'Q', 0, 'Water') / 1000.0
    if density_g_cm3 <= rho_sat:
        raise ValueError(f"Yoğunluk {density_g_cm3:.4f} g/cm3 doygun sıvı yoğunluğunun ({rho_sat:.4f}) altında")
    return PropsSI('T', 'D', density_g_cm3 * 1000.0, 'P', pressure_pa, 'Water')


def get_axial_density_profile(avg_value, gradient, num_axial_nodes=6, floor=0.0):
    n = int(num_axial_nodes)
    if n < 1:
        raise ValueError("num_axial_nodes >= 1 olmalı")
    return [max(avg_value * (1.0 + gradient * ((i + 0.5) / n - 0.5)), floor) for i in range(n)]


def fission_rate_density(power_density_w_cm3):
    return power_density_w_cm3 / ENERGY_PER_FISSION_J


def equilibrium_xenon_density(fission_rate, flux, sigma_a_xe_barn):
    n_per_cm3 = (GAMMA_I135 + GAMMA_XE135) * fission_rate / (LAMBDA_XE135 + sigma_a_xe_barn * 1.0e-24 * flux)
    return n_per_cm3 * 1.0e-24


def fresh_fuel_material(name='fuel'):
    m = openmc.Material(name=name)
    m.add_elements_from_formula('UO2', enrichment=ENRICHMENT_WO)
    m.set_density('g/cm3', FUEL_DENSITY)
    return m


def atom_densities(mat):
    return {n: float(v[1]) if isinstance(v, (tuple, list)) else float(v)
            for n, v in mat.get_nuclide_atom_densities().items()}


def _load_table():
    if 'table' not in _cache:
        if not os.path.exists(TABLE_PATH):
            raise FileNotFoundError(f"{TABLE_PATH} yok. Önce model_mutation.py'yi çalıştırın.")
        with open(TABLE_PATH) as f:
            _cache['table'] = json.load(f)
    return _cache['table']


def burnup_composition(burnup_gwd_t):
    if burnup_gwd_t <= 0.0:
        if 'fresh' not in _cache:
            _cache['fresh'] = atom_densities(fresh_fuel_material())
        return dict(_cache['fresh'])
    t = _load_table()
    grid, comps = t['burnup_gwd_t'], t['compositions']
    if burnup_gwd_t > grid[-1] + 1e-9:
        raise ValueError(f"burnup {burnup_gwd_t} tablo aralığının ({grid[-1]}) dışında")
    i = max(j for j, b in enumerate(grid) if b <= burnup_gwd_t + 1e-12)
    if i == len(grid) - 1:
        return dict(comps[i])
    w = (burnup_gwd_t - grid[i]) / (grid[i + 1] - grid[i])
    keys = set(comps[i]) | set(comps[i + 1])
    return {n: (1.0 - w) * comps[i].get(n, 0.0) + w * comps[i + 1].get(n, 0.0) for n in keys}


def make_fuel_material(name, burnup_gwd_t, xe135_atom_density):
    mat = openmc.Material(name=name)
    mat.set_density('sum')
    for nuc, dens in burnup_composition(burnup_gwd_t).items():
        if nuc not in EXCLUDED and dens > 0.0:
            mat.add_nuclide(nuc, dens)
    if xe135_atom_density > 0.0:
        mat.add_nuclide('Xe135', xe135_atom_density)
    return mat