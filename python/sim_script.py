import sys
import os
import argparse
import ucntrap_py

def run():
    # 0. Retrieve Dagger Mode
    parser = argparse.ArgumentParser(description="UCN Production Script")
    parser.add_argument('--mode', type=str, default='Fast',
                        choices=['Fast', 'Slow', 'Segmented'], help="Dagger detection mode")
    parser.add_argument('--out_dir', type=str, default='test',
                        help="Sub-folder under results/ (e.g., 'defect2e-3')")
    args, unknown = parser.parse_known_args()

    # 1. Initialize MPI
    ucntrap_py.init_mpi(sys.argv)

    # Use absolute path (will be ../ from current file)
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # 2. Obtain Slurm Array ID
    global_id = int(os.getenv('SLURM_ARRAY_TASK_ID', '0'))
    hold_times = [20.0, 50.0, 100.0, 200.0, 1550.0]

    ht_idx = global_id // 20
    batch_idx = global_id % 20
    current_ht = hold_times[ht_idx]

    # 3. Configure Simulation
    config = ucntrap_py.SimulationConfig()
    config.ntraj = 50000
    config.dt = 0.001
    config.seed = 42 + batch_idx
    config.field_model = "trap"
    config.hold_time = current_ht

    mode_map = {
        'Fast': ucntrap_py.DaggerMode.Fast,
        'Slow': ucntrap_py.DaggerMode.Slow,
        'Segmented': ucntrap_py.DaggerMode.Segmented
    }
    config.dagger_mode = mode_map[args.mode]

    # 4. Path fixing (absolute path)
    config.neutron_init_file = os.path.join(base_dir, "data/neutrons_init.out")
    config.x_trace_file = os.path.join(base_dir, "data/xvals.bin")
    config.y_trace_file = os.path.join(base_dir, "data/yvals.bin")
    config.z_trace_file = os.path.join(base_dir, "data/zvals.bin")

    # 5. Output structure
    output_dir = os.path.join(base_dir, "results", args.folder, f"HT_{int(current_ht)}_{args.mode}")
    os.makedirs(output_dir, exist_ok=True)
    config.output_prefix = os.path.join(output_dir, f"batch_{batch_idx:02d}")

    # 6. Execute
    print(f"Task {global_id}: HT={current_ht}s, Folder='{args.folder}', Mode={args.mode}")    with ucntrap_py.ostream_redirect():
        runner = ucntrap_py.Runner(config)
        runner.run()

if __name__ == "__main__":
    run()