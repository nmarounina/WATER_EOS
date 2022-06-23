from parameters import M_h2o
import EOS
import numpy as np
import math as m
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt




list_logp = np.linspace(3,10,100)
list_T = np.linspace(200, 3000, 100)
T_grid,logp_grid = np.meshgrid(list_T,list_logp)

rho_grid=np.zeros( (len(T_grid[:,0]), len(T_grid[0,:])) )

for i in range( 0,len(list_logp) ):
    for j in range( 0,len(list_T) ):

        point=EOS.DataPoint(T_grid[i,j], 10**logp_grid[i,j])
        rho_grid[i,j]=point.rho*M_h2o

        #print(" %f %e %f %s" %( point.T, point.p, point.rho, point.grid ))







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
