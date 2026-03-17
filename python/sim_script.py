import sys
import os
import ucntrap_py

def run_simulation():
    ucntrap_py.init_mpi(sys.argv)

    try:
        config = ucntrap_py.SimulationConfig()

        config.ntraj = 10000
        config.dt = 0.001
        config.seed = 42
        
        # --- 實驗參數 ---
        config.holding_time = 20.0
        config.cleaning_time = 50.0
        config.dagger_mode = ucntrap_py.DaggerMode.Fast

        # --- Directory Paths ---
        config.neutron_init_file = "../data/neutrons_init_test.out"
        config.x_trace_file = "../data/xvals.bin"
        config.y_trace_file = "../data/yvals.bin"
        config.z_trace_file = "../data/zvals.bin"
        
        # Output Prefix
        config.output_prefix = "../results/python_multi_rank/sim_run"

        # 3. Check folder existence
        output_dir = os.path.dirname(config.output_prefix)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        # 4. Initialize runner
        runner = ucntrap_py.Runner(config)

        # 5. Simulate
        print(f"Starting simulation with {config.ntraj} trajectories...")
        with ucntrap_py.ostream_redirect():
            status = runner.run()

        if status == 0:
            print(f"\nSimulation finished successfully. Results saved to: {config.output_prefix}.csv")
        else:
            print(f"\nSimulation exited with status: {status}")

    except RuntimeError as e:
        print(f"\nC++ Error encountered: {e}")
    except Exception as e:
        print(f"\nPython Error: {e}")

    return 0

if __name__ == "__main__":
    run_simulation()