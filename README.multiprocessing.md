# How to Run Parallel Fire Behavior Simulations

I've compared a variety of approaches to running parallel simulations using `pyretechnics.eulerian_level_set.spread_fire_with_phi_field`, and I've produced a modest collection of test scripts, which can be found in the toplevel `prof` folder of the repository.

Here is the high-level summary of the lessons I've learned from this exercise:

## tl;dr: The fastest and most memory efficient solution is:

- Use Python's native multiprocessing module.
- Use multiprocessing.set_start_method("spawn").
- Use multiprocessing.Pool(processes=$NUM_CORES) and pool.starmap(..., chunksize=max(1, $NUM_JOBS // $NUM_CORES)).
- Increase the jobs to cores ratio.
- Create all input arrays with dtype=numpy.float32.
- Use multiprocessing.shared_memory.SharedMemoryfor managing the input arrays.
- These techniques are all demonstrated together in prof/spread_fire_multiprocessing.py.

## In Detail

1. Use `dtype=numpy.float32` when creating arrays to be passed to `SpaceTimeCube` to avoid array copying.
2. A simple bash `for` loop, that runs a single fire per subprocess is the fastest solution, but duplicates all the input layers in memory (`for` + spread_fire.py).
3. The Python `multiprocessing` module performs around 250 ms slower than the `for` loop and uses the same amount of memory (spread_fire_multiprocessing_no_shared_mem.py).
4. Setting `multiprocessing.set_start_method("spawn")` and using a high chunksize within `multiprocessing.Pool()` provides a noticeable speed improvement and makes the code OS-agnostic.
5. The "Areal Throughput" metric reduces as the number of simultaneous cores (not jobs) increases although this is less apparent when the system's CPUs aren't saturated. This can be mitigated by increasing the ratio of jobs to cores and using a chunksize less than or equal to this ratio.
6. Using shared memory input arrays dramatically reduces memory usage (by around 5x) for around 750 ms of additional runtime in total (spread_fire_multiprocessing.py).
7. Using `multiprocessing.Pool`'s `initializer(*initargs)` feature provides no speed or memory improvements over not using it but creates resource leak errors (spread_fire_multiprocessing_init_worker.py).
8. The multithreading approach heavily underutilizes CPU cores (because of Python's GIL) and only uses about 5% less memory than multiprocessing with shared memory (spread_fire_multithreading.py).
9. The `dask` library uses ~7% more memory, ~20% more time, and reduces areal throughput by ~6%. This is even worse if a manual shared memory solution isn't implemented (spread_fire_dask.py).
10. The `ray.util.multiprocessing` module uses ~5 seconds of runtime and ~1-3+ GB of RAM more than Python's native `multiprocessing` module (spread_fire_ray_multiprocessing.py).
11. When using the `ray.put` function, make sure to run `ray.init(num_cpus=$NUM_CPUS)` early in the main thread in order to avoid a ~1 second runtime and ~3-5+ GB memory penalty over Python's native `multiprocessing.shared_memory` module (spread_fire_ray_multiprocessing_with_put.py).
12. There doesn't appear to be any performance difference between using `ray.remote` and `ray.util.multiprocessing` (spread_fire_ray_core_with_put.py).
13. When using a Ray cluster, your child processes are scheduled by default with low priority (e.g., `nice > 0`) and can be pre-empted by other higher priority processes.
