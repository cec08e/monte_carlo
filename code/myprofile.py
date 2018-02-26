import pstats, cProfile

import heisenberg3d_opt_cy

cProfile.runctx("heisenberg3d_opt_cy.sweep_k()", globals(), locals(), "Profile.prof")

s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()
