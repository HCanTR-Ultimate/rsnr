import os
import glob
import json
import itertools
import traceback
from datetime import datetime
import numpy as np
import openmc
from mutation import mutate_model, TEMPERATURE_SETTINGS

CONTROL_RODS_LEVELS = [0.0, 0.5, 1.0]
XE_AVG_LEVELS = [0.0, 3e-9, 6e-9, 9e-9]
XE_GRADIENT_LEVELS = [-0.3, 0.3]
COOLANT_DENSITY_LEVELS = [0.65, 0.75]
COOLANT_GRADIENT_LEVELS = [-0.05, 0.05]
BURNUP_LEVELS = [0.0, 10.0, 25.0, 40.0]

BATCHES = 340
INACTIVE = 40
PARTICLES = 10000

NUM_THREADS = "12"
NUM_NODES = 24
Z_HALF = 100.0
CORE_HALF = 75.25
SAVE_EVERY = 512


def get_combinations():
    combinations = []
    for banks in itertools.product(CONTROL_RODS_LEVELS, repeat=4):
        for xe_avg in XE_AVG_LEVELS:
            xe_gradients = XE_GRADIENT_LEVELS if xe_avg > 0.0 else [0.0]
            for xe_grad in xe_gradients:
                for rho in COOLANT_DENSITY_LEVELS:
                    for rho_grad in COOLANT_GRADIENT_LEVELS:
                        for burnup in BURNUP_LEVELS:
                            combinations.append((banks, xe_avg, xe_grad, rho, rho_grad, burnup))
    return combinations


def cleanup_output_files():
    for pattern in ['statepoint.*.h5', 'summary.h5', 'tallies.out', '*.xml']:
        for f in glob.glob(pattern):
            if f not in ['materials.xml', 'geometry.xml', 'plots.xml']:
                try:
                    os.remove(f)
                except OSError:
                    pass


def setup_settings_and_tallies():
    os.environ["OMP_NUM_THREADS"] = NUM_THREADS
    settings = openmc.Settings()
    settings.batches = BATCHES
    settings.inactive = INACTIVE
    settings.particles = PARTICLES
    settings.run_mode = 'eigenvalue'
    settings.temperature = TEMPERATURE_SETTINGS
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-CORE_HALF, -CORE_HALF, -Z_HALF),
                               (CORE_HALF, CORE_HALF, Z_HALF),
                               only_fissionable=True))

    mesh = openmc.CylindricalMesh(
        r_grid=[0.0, 100.0],
        phi_grid=[0.0, 2 * np.pi],
        z_grid=np.linspace(-Z_HALF, Z_HALF, NUM_NODES + 1)
    )
    tally = openmc.Tally(name='axial_tally')
    tally.filters = [openmc.MeshFilter(mesh)]
    tally.scores = ['fission', 'absorption']
    return settings, openmc.Tallies([tally])


def run_simulation(params):
    cr, xe_avg, xe_grad, rho_avg, rho_grad, burnup = params

    try:
        materials, geometry = mutate_model(cr, xe_avg, xe_grad, rho_avg, rho_grad, burnup)

        cleanup_output_files()
        materials.export_to_xml('materials.xml')
        geometry.export_to_xml('geometry.xml')

        settings, tallies = setup_settings_and_tallies()
        settings.export_to_xml('settings.xml')
        tallies.export_to_xml('tallies.xml')

        openmc.run(output=True)

        sp_file = glob.glob('statepoint.*.h5')[0]
        with openmc.StatePoint(sp_file) as sp:
            k_eff = float(sp.keff.nominal_value)
            k_std = float(sp.keff.std_dev)

            tally = sp.get_tally(name='axial_tally')
            mean = tally.mean.ravel()
            fiss_idx = tally.scores.index('fission')
            abs_idx = tally.scores.index('absorption')
            n_scores = len(tally.scores)

            fission_rates = mean[fiss_idx::n_scores][:NUM_NODES].tolist()
            absorption_rates = mean[abs_idx::n_scores][:NUM_NODES].tolist()

            total_fission = sum(fission_rates)
            if total_fission > 0:
                power_distribution = [f / total_fission for f in fission_rates]
            else:
                power_distribution = [0.0] * NUM_NODES

        cleanup_output_files()
        return {
            "status": "success",
            "k_eff": k_eff,
            "k_std": k_std,
            "fission_rates": fission_rates,
            "absorption_rates": absorption_rates,
            "power_distribution": power_distribution
        }
    except Exception as e:
        cleanup_output_files()
        traceback.print_exc()
        return {"status": "error", "message": str(e)}


def main():
    combinations = get_combinations()
    total_comb = len(combinations)
    print(f"Toplam Kombinasyon Sayısı: {total_comb}")

    range_input = input("Hesaplanacak aralık (default: -1, örn: 100-120): ").strip()
    if not range_input or range_input == "-1":
        indices = range(total_comb)
    else:
        try:
            start_idx, end_idx = map(int, range_input.split('-'))
            indices = range(max(0, start_idx), min(total_comb, end_idx + 1))
        except ValueError:
            print("Geçersiz format, tümü hesaplanıyor.")
            indices = range(total_comb)

    target_count = len(indices)

    proceed = input(f"Seçilen {target_count} kombinasyon için işlem başlatılacak. Devam edilsin mi? (y/n): ").strip().lower()
    if proceed != 'y':
        print("İşlem iptal edildi.")
        return

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    session_dir = os.path.join("save", f"run_{timestamp}_{target_count}comb")
    os.makedirs(session_dir, exist_ok=True)

    results = []
    for i, idx in enumerate(indices, 1):
        params = combinations[idx]
        res = run_simulation(params)

        results.append({
            "index": idx,
            "parameters": {
                "control_rods": list(params[0]),
                "xe_avg": params[1],
                "xe_gradient": params[2],
                "coolant_density_avg": params[3],
                "coolant_gradient": params[4],
                "burnup": params[5]
            },
            "result": res
        })

        print(f"[{i}/{target_count}] Kombinasyon {idx} tamamlandı. Durum: {res.get('status')}")

        if i % SAVE_EVERY == 0 or i == target_count:
            batch_filename = os.path.join(session_dir, f"results_batch_{i}.json")
            with open(batch_filename, 'w', encoding='utf-8') as f:
                json.dump({
                    "timestamp": timestamp,
                    "total_combinations": target_count,
                    "completed_count": i,
                    "data": results
                }, f, indent=4)

    final_filename = os.path.join(session_dir, "results_final.json")
    with open(final_filename, 'w', encoding='utf-8') as f:
        json.dump({
            "timestamp": timestamp,
            "total_combinations": target_count,
            "data": results
        }, f, indent=4)

    print(f"Tüm hesaplamalar tamamlandı. Kayıt konumu: {session_dir}")


if __name__ == "__main__":
    main()