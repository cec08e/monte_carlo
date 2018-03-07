from heisenberg3d_opt_cy import sweep_k
import time

start = time.time()
sweep_k()
end = time.time()
print("Execution time (Cython):", end-start)
