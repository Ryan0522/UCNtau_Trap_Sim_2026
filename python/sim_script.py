import sys
import os
import argparse
sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), 'python'))
import ucntrap_py
import resource

def run():
    # 0. Retrieve Dagger Mode
    parser = argparse.ArgumentParser(description="UCN Production Script")
    parser.add_argument('--mode', type=str, default='Fast',
                        choices=['Fast', 'Slow', 'Segmented'], help="Dagger detection mode")
    parser.add_argument('--defect', type=float, default=0,
                        help="Magnetic Field Defect Parameter")
    parser.add_argument('--out_dir', type=str, default='test',
                        help="Sub-folder under results/ (e.g., 'defect2e-3')")
    parser.add_argument('--ntraj', type=int, default=100, 
                        help="# of simulations")
    args, unknown = parser.parse_known_args()

    # 1. Initialize MPI
    ucntrap_py.init_mpi(sys.argv)

    try:
        # Use absolute path (will be ../ from current file)
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        # 2. Obtain Slurm Array ID
        global_id = int(os.getenv('SLURM_ARRAY_TASK_ID', '0'))
        total_ntraj_per_task = args.ntraj
        hold_times = [20.0, 50.0, 100.0, 200.0, 1550.0]

        rank = 0
        try:
            rank = untrap_py.get_rank()
        except:
            rank = int(os.getenv('OMPI_COMM_WORLD_RANK', '0'))

        ht_idx = global_id // 4
        batch_idx = global_id % 4
        current_ht = hold_times[ht_idx]

        # 3. Configure Simulation
        config = ucntrap_py.SimulationConfig()
        config.ntraj = args.ntraj
        config.dt = 0.001
        config.defect = args.defect
        config.seed = 42 + (global_id * 1000) + rank
        config.holding_time = current_ht
        config.dip_end_times = [current_ht, current_ht + 250.0]
        config.array_offset = global_id * total_ntraj_per_task
        config.heat_mult = 0.0 # For energy conservation test

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
        output_dir = os.path.join(base_dir, "results", args.out_dir, f"HT_{int(current_ht)}_{args.mode}")
        os.makedirs(output_dir, exist_ok=True)
        config.output_prefix = os.path.join(output_dir, f"batch_{batch_idx:02d}")

        # 6. Execute
        print(f"Task {global_id}: HT={current_ht}s, Dir='{args.out_dir}', Mode={args.mode}, Defect={args.defect}")
        with ucntrap_py.ostream_redirect():
            runner = ucntrap_py.Runner(config)
            runner.run()
    
    except Exception as e:
        print(f"Error in Rank Task {global_id}: {e}")
    
    finally:
        if hasattr(ucntrap_py, 'finalize_mpi'):
            ucntrap_py.finalize_mpi()
        else:
            print("Warning: finalize_mpi not found. Please recompile C++ wrapper.")
        
        try:
            peak_mem_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            peak_mem_mb = peak_mem_kb / 1024.0
            print(f"\n--- Resource Usage Task {global_id} Rank {rank} ---")
            print(f"Peak Memory Usage: {peak_mem_mb:.2f} MB")
        except Exception as mem_e:
            print(f"Could not retrieve memory usage: {mem_e}")


if __name__ == "__main__":
    run()