from heisenberg3d_opt_cy import sweep_k, plot_M_v_B, Heisenberg3D
import time

start = time.time()
#sweep_k()
plot_M_v_B()
#lat = Heisenberg3D(10, 10, init_T=5, J_inter=.1)
#lat.mag_v_temp(init_temp=5, final_temp = .01, temp_step =.05, eq_time=7000, cor_time=1000)
end = time.time()
print("Execution time (Cython):", end-start)
