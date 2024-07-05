#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
ny = 8
nx = 16
poly_degree = 5
# Re_tau = 5200.0
# Re_tau = 395.0
Re_tau = 5186.0
# ========================================
delta = 1.0
L_inf = delta
Lx = 2.0*np.pi*delta
Ly = 2.0*delta
Lz = 3.141592653589793238*delta

y_distance = Ly/np.float64(ny)
dy_effective = y_distance/np.float64(poly_degree)
dx_effective = Lx/np.float64(nx*poly_degree)
print(y_distance)
print(dy_effective)
print(dx_effective)

mu_ref = 1.716e-5
T_ref = 273.15
density = 1.225 # kg/m3

frere_dy_plus_effective = 250.0
frere_u_friction = mu_ref*Re_tau/density/delta
dy_plus_calc = dy_effective*frere_u_friction*density/mu_ref
dx_plus_calc = dx_effective*frere_u_friction*density/mu_ref
print(dy_plus_calc)
print(dx_plus_calc)
y_plus_calc = y_distance*frere_u_friction*density/mu_ref
print(y_plus_calc)

tau_wall_from_Re_tau = pow(mu_ref*Re_tau/delta,2.0)/density

Re_bulk = np.power(0.073, -4.0/7.0)*np.power(2.0, 5.0/7.0)*np.power(Re_tau, 8.0/7.0)
print("Re bulk:")
print(Re_bulk)
velocity_bulk = mu_ref*Re_bulk/(density*delta)
latest_tau_w_value_from_simulation = 2.6483112826572973e-04#1.4602213022727240e-04
tau_wall_from_code = latest_tau_w_value_from_simulation*(density*velocity_bulk*velocity_bulk)
print("tau_wall_from_Re_tau: ")
print(tau_wall_from_Re_tau)
print("tau_wall_from_code: ")
print(tau_wall_from_code)

Cf = 0.073*np.power(2.0*Re_bulk, -1.0/4.0)
print("Re bulk:")
print(Re_bulk)
print("time two:")
print(2.0*Re_bulk)
print("Cf value w/o modifications")
print("%1.3e" % Cf)
Cf_per_area = Cf/(Lx*Lz)
val_to_print = 2.0*Cf_per_area
print("Cf value w/ modifications")
print("%1.3e" % val_to_print)
# Cf = 0.073*np.power(14130.0, -1.0/4.0)
# print("Lodato's cf:")
# print(Cf)

# Iter: 195    Time: 0.0733571    Cf: 0.000690282    Ub: 0.961797    BulkMassFlow: 0.961797