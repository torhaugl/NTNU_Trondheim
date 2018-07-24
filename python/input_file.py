# SHORT
final_time = 14 * 60.0  # minutes
conc_bulk = 0.2
dt = 0.5 / 60  # min = 0.5 sec [5 sec]

# LONG
N = int(final_time / dt)  # Number of steps

vortex_length = vl = 17.0  # Length of lattice (micro m), square
Lx, Ly, Lz = 3 * vl, 3 * vl, 10 * vl  # Length of Biofilm

max_mass_cell = 14700
avg_mass_cell = 410
density_cell = 290
density_eps = 290
max_mass_eps = density_eps / density_cell * max_mass_cell

max_volume_fraction = 0.52
diffusion_subst = 40680
diffusion_qsm = diffusion_qsi = 33300
max_particles = density_cell * max_volume_fraction * vl ** 3 / max_mass_cell

half_saturation = Ks = 2.34e-3
max_substrate_uptake_rate = Vmax = 0.046
max_yield = 0.444
maintenance_rate = 6e-4

Zed = 0  # Down-regulated EPS production
Zeu = 1e-3  # Up-regulated EPS production
Zqd = 8.3
Zqu = 1230
transfer_coefficient = 1e-3

alpha = 1.33
beta = 10
gamma = 0.001
Kq = 10  # Positive feedback (QSM production with QSM concentration)
