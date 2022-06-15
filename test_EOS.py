import numpy as np

from parameters import M_h2o, Rig
import IAPWS95
from scipy import optimize
import numpy as np
import math as m
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt



def get_a_datatpint_from_EOS(T, p):
    point=IAPWS95.DataPoint( T, p )
    print( " %s %f %s" %("Density =",  point.rho*M_h2o, "kg.m-3" ))
    print(" %s %f %s" %("Entropy =" , point.s, "J.mol-1.K-1"))
    print(" %s %f %s" %("Heat capacity =" , point.Cp, "J.mol-1.K-1"))
    print(" %s %f %s" %("Internal energy =" , point.u, "J.mol-1"))
    return 0



list_logp = np.linspace(3,9,100)
list_T = np.linspace(300, 1400, 50)
T_grid,logp_grid = np.meshgrid(list_T,list_logp)

rho_grid=np.zeros( (len(T_grid[:,0]), len(T_grid[0,:])) )

for i in range( 0,len(list_logp) ):
    for j in range( 0,len(list_T) ):

        point=IAPWS95.DataPoint(T_grid[i,j], 10**logp_grid[i,j])
        rho_grid[i,j]=point.rho*M_h2o


fig=plt.figure()

ax = plt.axes(projection='3d')

ax.scatter3D(647.096, m.log10(22.064e6), 322.,
             c="red",
             marker="*",
             s=80
             )

ax.plot_surface(T_grid, logp_grid, rho_grid,
                cmap="viridis",
                alpha=0.8)
ax.set_title("Density surface and the pure water critical point")
ax.set_xlabel("Temperature, K")
ax.set_ylabel("log10(Pressure, Pa)")
ax.set_zlabel("Density, kg.m$^{-3}$")
plt.show()


#Benchmark, see table 6.6 of Wagner and Pruss (2002)
#T=500 K and rho=838.025 kg.m-3 :
# def func_search_for_p(p):
#
#     point=IAPWS95.DataPoint( 500., p )
#     #print( "%e %f" %(p, point.rho*M_h2o) )
#
#     return (point.rho*M_h2o - 838.025)/838.025
#
# pressure=optimize.bisect(func_search_for_p, 1e7, 2.63e8)
#
#
#
# point_benchmark=IAPWS95.DataPoint( 500., pressure)
#
# print("\n")
# print( "Benchmark:")
# print( "T=", point_benchmark.T, "rho=", point_benchmark.rho*M_h2o)
# print("\n")
# print("Ideal-gas part and its derivatives:")
# print( point_benchmark.phi0 )
# print( point_benchmark.dphi0t )
# print( point_benchmark.dphi0tt )
#
# print("\n")
# print("Residual part and its derivatives:")
# print( point_benchmark.phir )
# print( point_benchmark.dphird )
# print( point_benchmark.dphirt )
# print( point_benchmark.dphirdd )
# print( point_benchmark.dphirtt )
# print( point_benchmark.dphirdt )
