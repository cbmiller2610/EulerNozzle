from ctypes import *
import numpy as np
import timeit
#Import C dlls
heat_solver_dll = CDLL("./heat_solver.so")
heat_solver_dll.argtypes = [c_double, c_int, c_double, c_int, c_double, c_double, c_double, c_double]
#set up grid
xmax = 1
xnodes = 301
dx = xmax/(xnodes-1)
#t=900
#tmax = t*1.60E-07/.15 # thermal diffusivity of water near 350K, and characteristic length 0.15m (6in)
#print("Simulation Time (s)")
#print(tmax)
tmax = 2
tnodes = 180001
dt = tmax/(tnodes-1)
#set up physical values
T0 = 350 # initial temp in Kelvin, a HOT cup of coffee (170deg F).
Tinf = 295 # room temperature air for convection
U0 = (T0-Tinf)/(T0-Tinf)
r = dt/(dx**2)
print("r value")
print(r)
Nu = 40 #educated guess for appropriate Nusselt Number with the delta between T0 and Tinf
#Create main results array
U_size = tnodes*xnodes
U_type = c_double * U_size
U = U_type()
U_ptr = byref(U)
U_np_1D = np.ctypeslib.as_array(U)
U_np_2D = np.reshape(U_np_1D, (tnodes,xnodes))
U_np_2D[0,:] = U0
#heat flux from electric heating
#C = -1.25
C = -2.5
status = heat_solver_dll.heat_solver(U_ptr, c_int(xnodes), c_double(dx), c_int(tnodes), c_double(dt),
c_double(Nu), c_double(r), c_double(C))
with open('301_spatial_180000tnodes_case3_2xC.dat','wb') as f:
    np.save(f, U_np_2D)
    #np.savetxt("U_results.csv",U_np_2D, "%10.10f", ',')
