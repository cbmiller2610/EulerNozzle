import numpy as np
import pdb
#from ctypes import *
import timeit

#Functions
def grid_transform(xi, eta, m1, m2, nozzle_xloc, nozzle_ymax, throat_xloc, throat_ymax):

    #Region 1: converging
    if xi <= 1.5:
        y = eta * (m1 * (xi - inlet_xloc) + inlet_ymax)
    #Region 2: diverging
    elif (xi > 1.5):
        y = eta*(m2 * (xi - throat_xloc) + throat_ymax)
    else:
        print("OUT OF RANGE")
    return y

def initial_cond(xi):
    if (xi >= 0) & (xi < 0.5):
        rho = 1.0
        T = 1.0
    elif (xi >= 0.5) & (xi < 1.5):
        rho = 1.0 - 0.366 * (xi - 0.5)
        T = 1.0 - 0.167 * (xi - 0.5)
    elif (xi >= 1.5) & (xi < 2.1):
        rho = 0.634 - 0.702 * (xi - 1.5)
        T = 0.833 - 0.4908 * (xi - 1.5)
    elif (xi >= 2.1) & (xi <= 3.0):
        rho = 0.5892 + 0.10228 * (xi - 2.1)
        T = 0.93968 + 0.0622 * (xi - 2.1) 
    else:
        print("OUT OF RANGE")
    return rho, T

def U_to_prim(U, xi_div, eta_div):
    #returns rho, u, v, and p primitive variable arrays
    U1 = U[:, :, 0]
    U2 = U[:, :, 1]
    U3 = U[:, :, 2]
    U4 = U[:, :, 3]
    U_step = np.empty((xi_div, eta_div, 5))
    g = 1.4
    U_step[:, :, 0] = U1 #rho
    U_step[:, :, 1] = U2/U1 #u
    U_step[:, :, 2] = U3/U1 #v
    U_step[:, :, 3] = (g - 1) * (U4 - (g / 2) * ((U2**2 + U3**2) / U1)) #pressure
    U_step[:, :, 4] = (g - 1) * (U4 / U1 - (g / 2) * ((U2**2 + U3**2) / U1**2)) #temperature
    return U_step

def U_to_E(U, xi_div, eta_div):
    U1 = U[:, :, 0]
    U2 = U[:, :, 1]
    U3 = U[:, :, 2]
    U4 = U[:, :, 3]
    E_step = np.empty((xi_div, eta_div, 4))
    g = 1.4
    E_step[:, :, 0] = U2
    E_step[:, :, 1] = U2**2 / U1 + (1 - 1 / g) * (U4 - (g / 2) * ((U2**2 + U3**2) / U1))
    E_step[:, :, 2] = (U2 * U3) / U1
    E_step[:, :, 3] = (g * U2 * U4) / U1 - ((g * (g - 1)) / 2) * ((U2**3 + U2 * U3**2)/U1)
    return E_step

def U_to_F(U, xi_div, eta_div):
    U1 = U[:, :, 0]
    U2 = U[:, :, 1]
    U3 = U[:, :, 2]
    U4 = U[:, :, 3]
    F_step = np.empty((xi_div, eta_div, 4))
    g = 1.4
    F_step[:, :, 0] = U3
    F_step[:, :, 1] = (U2 *U3) / U1
    F_step[:, :, 2] = U3**2 / U1 + (1 - 1 / g) * (U4 - (g / 2) * ((U2**2 + U3**2) / U1))
    F_step[:, :, 3] = (g * U3 * U4) / U1 - ((g * (g - 1)) / 2) * ((U2**2 * U3 + U3**3) / U1)
    return F_step

def Ughost_to_Eghost(U, xi_div):
    U1 = U[:, :, 0]
    U2 = U[:, :, 1]
    U3 = U[:, :, 2]
    U4 = U[:, :, 3]
    E_step = np.empty((xi_div-2, 2, 4))
    g = 1.4
    E_step[:, :, 0] = U2
    E_step[:, :, 1] = U2**2 / U1 + (1 - 1 / g) * (U4 - (g / 2) * ((U2**2 + U3**2) / U1))
    E_step[:, :, 2] = (U2 * U3) / U1
    E_step[:, :, 3] = (g * U2 * U4) / U1 - ((g * (g - 1)) / 2) * ((U2**3 + U2 * U3**2)/U1)
    return E_step

def Ughost_to_Fghost(U, xi_div):
    U1 = U[:, :, 0]
    U2 = U[:, :, 1]
    U3 = U[:, :, 2]
    U4 = U[:, :, 3]
    F_step = np.empty((xi_div-2, 2, 4))
    g = 1.4
    F_step[:, :, 0] = U3
    F_step[:, :, 1] = (U2 *U3) / U1
    F_step[:, :, 2] = U3**2 / U1 + (1 - 1 / g) * (U4 - (g / 2) * ((U2**2 + U3**2) / U1))
    F_step[:, :, 3] = (g * U3 * U4) / U1 - ((g * (g - 1)) / 2) * ((U2**2 * U3 + U3**3) / U1)
    return F_step

#Import C dlls
#solver_dll = CDLL("./euler_solver.so")
#solver_dll.argtypes = [c_double, c_double, c_double, c_int, c_int, c_int]

#Gridding Inputs
xi_max = 3
eta_max = 1
xi_div = 201
eta_div = 101
dxi = xi_max / (xi_div - 1)
deta = eta_max / (eta_div - 1)
dt = 0.0008
tsteps = 1000

#Grid Generation

xi_vec = np.linspace(0, xi_max, xi_div)
eta_vec = np.linspace(0, eta_max, eta_div)

xi_grid,eta_grid = np.meshgrid(xi_vec, eta_vec, indexing='ij')
x_grid = xi_grid

grid_trans_func = np.vectorize(grid_transform)
throat_area = 1
throat_diam = 2 * np.sqrt(throat_area / np.pi)
throat_ymax = throat_diam / 2
inlet_area = 5.95
nozzle_area= 5.95
inlet_diam = 2 * np.sqrt(inlet_area / np.pi)
inlet_ymax = inlet_diam / 2
nozzle_diam = 2 * np.sqrt(nozzle_area / np.pi)
nozzle_ymax = nozzle_diam / 2

throat_xloc = 1.5
inlet_xloc = 0
nozzle_xloc = 3

m1 = (throat_ymax - inlet_ymax) / (throat_xloc - inlet_xloc)
m2 = (nozzle_ymax - throat_ymax) / (nozzle_xloc - throat_xloc)
y_grid = grid_trans_func(xi_grid, eta_grid, m1, m2, nozzle_xloc, nozzle_ymax, throat_xloc, throat_ymax)

#Initial Conditions

IC_func = np.vectorize(initial_cond)
rho, T = IC_func(xi_grid)
u = 0.59 / rho
v = np.zeros((xi_div, eta_div))
g = 1.4

#Creation of the Tensors of Conserved Variables and Fluxes as C Arrays
#tensize = xi_div * eta_div * 4
#tentype = c_double * tensize
#U = tentype()
#U_ptr = byref(U)
#E = tentype()
#E_ptr = byref(E)
#F = tentype()
#F_ptr = byref(F)
#U_np_1D = np.ctypeslib.as_array(U)
#U_np_3D = np.reshape(U_np_1D, (xi_div, eta_div, 4))
#E_np_1D = np.cytypeslib.as_array(E)
#E_np_3D = np.reshape(E_np_1D, (xi_div, eta_div, 4))
#F_np_1D = np.cytypeslib.as_array(F)
#F_np_3D = np.reshape(F_np_1D, (xi_div, eta_div, 4))
P_stor = np.empty((xi_div, eta_div, 5, tsteps+1)) #store rho, u, v, p, and T values
U_np_3D = np.empty((xi_div, eta_div, 4))
E_np_3D = np.empty((xi_div, eta_div, 4))
F_np_3D = np.empty((xi_div, eta_div, 4))

U1 = rho
U2 = rho * u
U3 = rho * v
U4 = rho * (T / (g - 1) + (g / 2) * (u**2 + v**2))

U_np_3D[:, :, 0] = U1
U_np_3D[:, :, 1] = U2
U_np_3D[:, :, 2] = U3
U_np_3D[:, :, 3] = U4

E_np_3D[:, :, 0] = U2
E_np_3D[:, :, 1] = U2**2 / U1 + (1 - 1 / g) * (U4 - (g / 2) * ((U2**2 + U3**2) / U1))
E_np_3D[:, :, 2] = (U2 * U3) / U1
E_np_3D[:, :, 3] = (g * U2 * U4) / U1 - ((g * (g - 1)) / 2) * ((U2**3 + U2 * U3**2) / U1)

F_np_3D[:, :, 0] = U3
F_np_3D[:, :, 1] = (U2 * U3) / U1
F_np_3D[:, :, 2] = U3**2 / U1 + (1 - 1 / g) * (U4 - (g / 2) * ((U2**2 + U3**2) / U1))
F_np_3D[:, :, 3] = (g * U3 * U4) / U1 - ((g * (g - 1)) / 2) * ((U2**2 * U3 + U3**3) / U1)

P_stor[:, :, :, 0] = U_to_prim(U_np_3D, xi_div, eta_div)
Ughosts = np.empty((xi_div-2, 2, 4))
Eghosts = np.empty((xi_div-2, 2, 4))
Fghosts = np.empty((xi_div-2, 2, 4))

for i in range(1,tsteps):
    #Update boundary nodes from last timestep
    #left boundary, subsonic inlet (2 float, 2 prescribed)
    T = 1
    rho = 1
    U_np_3D[0, :, 0] = rho
    U_np_3D[0, :, 1:3] = 2 * U_np_3D[1, :, 1:3] - U_np_3D[2, :, 1:3]
    u0 = U_np_3D[0, :, 1] / U_np_3D[0, :, 0]
    v0 = U_np_3D[0, :, 2] / U_np_3D[0, :, 0]
    #U_prim = U_to_prim(U_np_3D, xi_div, eta_div)
    U_np_3D[0, :, 3] = rho * (T / (g - 1) + (g / 2) * (u0**2 + v0**2))

    #right boundary, supersonic outlet, all quantities float
    U_np_3D[-1, :, :] = 2 * U_np_3D[-2, :, :] - U_np_3D[-3, :, :]
    U_init = U_np_3D
    #Wall and symmetry boundaries (update ghosts)
    #Lower Boundary
    Ughosts[:, 0, :] = U_np_3D[1:-1, 0, :]
    Ughosts[:, 0, 2] = - U_np_3D[1:-1, 0, 2]
    #Upper Boundary
    Ughosts[:, 1, :] = U_np_3D[1:-1, -1, :]
    Ughosts[:, 1, 2] = - U_np_3D[1:-1, -1, 2]
    
    #Update Flux Terms
    E_np_3D = U_to_E(U_np_3D, xi_div, eta_div)
    F_np_3D = U_to_F(U_np_3D, xi_div, eta_div)
    #Eghosts = Ughost_to_Eghost(Ughosts, xi_div)
    Fghosts = Ughost_to_Fghost(Ughosts, xi_div)

    #Predictor Step (Upper boundary uses ghost)
   #Update Upper Before Internal Flow
    U_np_3D[1:-1, -1, :]  = (U_np_3D[1:-1, -1, :] - (dt / dxi) * (E_np_3D[2:xi_div, -1, :] - E_np_3D[1:-1, -1, :])
                            - (dt / deta) * (Fghosts[:, 1, :] - F_np_3D[1:-1, -1, :]))
    #Inner
    U_np_3D[1:-1, 0:-1, :] = (U_np_3D[1:-1, 0:-1, :]
                            - (dt / dxi) * (E_np_3D[2:xi_div, 0:-1, :] - E_np_3D[1:-1, 0:-1, :])
                            - (dt / deta) * (F_np_3D[1:-1, 1:eta_div, :] - F_np_3D[1:-1, 0:-1, :]))

    #Corrector Step (Lower boundary uses ghost)
    #Update Flux Terms to Next Time Level
    E_np_3D = U_to_E(U_np_3D, xi_div, eta_div)
    F_np_3D = U_to_F(U_np_3D, xi_div, eta_div)
    #Update Lower Before Internal Flow
    U_np_3D[1:-1, 0, :]  = 0.5 * (U_np_3D[1:-1, 0, :] + U_init[1:-1, 0, :]
                            - (dt / dxi) * (E_np_3D[1:-1, 0, :] - E_np_3D[0:-2, 0, :])
                            - (dt / deta) * (F_np_3D[1:-1, 0, :] - Fghosts[:, 0, :]))
    #Inner
    U_np_3D[1:-1, 1:eta_div, :] = 0.5 * (U_np_3D[1:-1, 1:eta_div, :]  + U_init[1:-1, 1:eta_div, :]
                                         - (dt / dxi) * (E_np_3D[1:-1, 1:eta_div, :] - E_np_3D[0:-2, 1:eta_div, :])
                                         - (dt / deta) * (F_np_3D[1:-1, 1:eta_div, :] - F_np_3D[1:-1, 0:-1, :]))
    
    #Transform Back to Physical Coordinates
    xi_max = 3
    eta_max = 1
    xi_vec = np.linspace(0,xi_max,xi_div)
    eta_vec = np.linspace(0,eta_max,eta_div)
    xi_grid,eta_grid = np.meshgrid(xi_vec,eta_vec,indexing='ij')


    P_stor[:, :, :, i] = U_to_prim(U_np_3D, xi_div, eta_div)

with open('TEST.dat', 'wb') as f: 
    np.save(f, P_stor)
