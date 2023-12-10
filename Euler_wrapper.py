import numpy as np
import pdb
from ctypes import *
import timeit

#Functions
def grid_transform(xi, eta, m1, m2, nozzle_xloc,nozzle_ymax,throat_xloc,throat_ymax):

    #Region 1: converging
    if xi<=1.5:
        y=eta*(m1*(xi-inlet_xloc)+inlet_ymax)
    #Region 2: diverging
    elif (xi > 1.5):
        y=eta*(m2*(xi-throat_xloc)+throat_ymax)
   else:
        print("WRONG VALUE")
    return y

#Import C dlls
solver_dll = CDLL("./euler_solver.so")
solver_dll.argtypes = [c_double, c_int, c_double, c_int, c_double, c_double, c_double, c_double]

#Gridding Inputs
xi_max = 3
eta_max = 1
xi_div = 100
eta_div = 50
dt=0.0015

#Grid Generation

xi_vec = np.linspace(0,xi_max,xi_div)
eta_vec = np.linspace(0,eta_max,eta_div)

xi_grid,eta_grid = np.meshgrid(xi_vec,eta_vec,indexing='ij')
x_grid = xi_grid

grid_trans_func = np.vectorize(grid_transform)
throat_area = 1
throat_diam = 2*np.sqrt(throat_area/np.pi)
throat_ymax = throat_diam/2
inlet_area = 5.95
nozzle_area= 5.95
inlet_diam = 2*np.sqrt(inlet_area/np.pi)
inlet_ymax = inlet_diam/2
nozzle_diam = 2*np.sqrt(nozzle_area/np.pi)
nozzle_ymax = nozzle_diam/2

throat_xloc = 1.5
inlet_xloc = 0
nozzle_xloc = 3

m1 = (throat_ymax - inlet_ymax)/(throat_xloc - inlet_xloc)
m2 = (nozzle_ymax - throat_ymax)/(nozzle_xloc - throat_xloc)
y_grid = grid_trans_func(xi_grid,eta_grid, m1, m2, nozzle_xloc,nozzle_ymax,throat_xloc,throat_ymax)

