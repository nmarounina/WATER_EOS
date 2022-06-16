from parameters import M_h2o
import math as m
import numpy as np
from scipy import interpolate
from scipy import constants as cst


def create_grids_for_interpolation_on_Mazevet_data():
    Mazevet_raw = np.loadtxt("./dataIN/eosMazevet.dat")  ### this file should be a SQUARE GRID IN T-rho

    # Columns in the Mazevet file:
    # rho[g/cc], T [K], P [GPa], P/(n_i kT), F/(N_i kT), U/(N_i kT), C_V/(N_i k), chi_T, chi_r, S/(N_i k)
    s_offset = -54.280936454849495  # J.mol-1.K-1

    # Grid 1 = square grid in T, rho
    RHO_raw = Mazevet_raw[:, 0]
    P_raw = Mazevet_raw[:, 2]
    T_raw = Mazevet_raw[:, 1]
    S_raw = Mazevet_raw[:, 9]
    U_raw = Mazevet_raw[:, 5]


    # let's do the interpolation on the p, T square grid
    N_interpol_grid = 20  # size of the new square grid
    logp = np.linspace(m.log10(1e-1),
                       m.log10(1e4),
                       N_interpol_grid)  # GPa
    p_interpol_grid = 10 ** logp
    T_interpol_grid = np.linspace(m.log10(np.min(T_raw + 5)),
                                  m.log10(np.max(T_raw - 5)),
                                  N_interpol_grid)
    T_interpol_grid = 10 ** T_interpol_grid  # K

    T_meshgrid, p_meshgrid = np.meshgrid(T_interpol_grid, p_interpol_grid)

    S_forRBS_Tp = interpolate.griddata((T_raw, P_raw), S_raw, (T_meshgrid, p_meshgrid), method="linear")
    U_forRBS_Tp = interpolate.griddata((T_raw, P_raw), U_raw, (T_meshgrid, p_meshgrid), method="linear")
    RHO_forRBS_Tp = interpolate.griddata((T_raw, P_raw), RHO_raw, (T_meshgrid, p_meshgrid), method="linear")

    # unit correction::
    for i in range(0, N_interpol_grid):
        for j in range(0, N_interpol_grid):
            S_forRBS_Tp[i, j] = S_forRBS_Tp[i, j] * cst.k * 3. * cst.N_A + s_offset
            U_forRBS_Tp[i,j]  = U_forRBS_Tp[i, j] * cst.k * 3. * cst.N_A * T_meshgrid[i,j]
            RHO_forRBS_Tp[i,j]= RHO_forRBS_Tp[i,j] * 1000. / M_h2o


    ######################### INSERT THOSE FUNCTIONS IN A FUTURE RELEVANT SPOT
    # RBS require 2 1D vectors for x and y, and a "meshgridded" 2D vector for Z
    s_mazevet_Tp = interpolate.RectBivariateSpline(p_interpol_grid * 1e9,
                                                   T_interpol_grid,
                                                   S_forRBS_Tp)

    # T in K, p in Pa, gives U in J.mol-1
    u_mazevet_Tp = interpolate.RectBivariateSpline(p_interpol_grid * 1e9,
                                                   T_interpol_grid,
                                                   U_forRBS_Tp )

    # T in K, p in Pa, gives RHO in mol.m-3
    rho_mazevet_Tp = interpolate.RectBivariateSpline(p_interpol_grid * 1e9,
                                                     T_interpol_grid,
                                                     RHO_forRBS_Tp)

    return p_interpol_grid * 1e9, T_interpol_grid, S_forRBS_Tp, U_forRBS_Tp, RHO_forRBS_Tp


def create_grids_for_interpolation_on_HpLT_grid():
    HPLT_raw = np.loadtxt("./dataIN/dataHPices.dat")

    N_lines_raw = len(HPLT_raw)
    N_square_raw_grid =int(N_lines_raw ** 0.5)
    alpha_interpol_grid = np.zeros(( N_square_raw_grid, N_square_raw_grid))
    s_interpol_grid = np.zeros((N_square_raw_grid, N_square_raw_grid))
    Cp_interpol_grid = np.zeros((N_square_raw_grid, N_square_raw_grid))
    rho_interpol_grid = np.zeros((N_square_raw_grid, N_square_raw_grid))

    T_raw = HPLT_raw[:, 0]
    p_raw = HPLT_raw[:, 1]
    rhoHPI = HPLT_raw[:, 2]
    alpha_raw = HPLT_raw[:, 3]
    Cp_raw = HPLT_raw[:, 4]
    s_raw = HPLT_raw[:, 5]

    T_interpol_grid = np.linspace(np.min(T_raw),
                                  np.max(T_raw),
                                  N_square_raw_grid)
    p_interpol_grid = np.linspace(m.log10(np.min(p_raw)),
                                  m.log10(np.max(p_raw)),
                                  N_square_raw_grid)
    p_interpol_grid = 10 ** p_interpol_grid


    cpt = 0
    for i in range(0, N_square_raw_grid):
        for j in range(0, N_square_raw_grid):
            alpha_interpol_grid[i, j] = alpha_raw[cpt]
            s_interpol_grid[i, j] = s_raw[cpt]
            Cp_interpol_grid[i, j] = Cp_raw[cpt]
            rho_interpol_grid[i, j] = rho_raw[cpt]

            cpt += 1

    ################ INSERT IN AN APPROPRIATE PLACE LATER
    s_HPI_Tp = interpolate.RectBivariateSpline(p_interpol_grid,
                                               T_interpol_grid,
                                               s_interpol_grid)
    alpha_HPI_Tp = interpolate.RectBivariateSpline(p_interpol_grid,
                                                   T_interpol_grid,
                                                   alpha_interpol_grid)
    rho_HPI_Tp = interpolate.RectBivariateSpline(p_interpol_grid,
                                                 T_interpol_grid,
                                                 rho_interpol_grid)
    Cp_HPI_Tp = interpolate.RectBivariateSpline(p_interpol_grid,
                                                T_interpol_grid,
                                                Cp_interpol_grid)

    return p_interpol_grid,T_interpol_grid,rho_interpol_grid,s_interpol_grid,Cp_interpol_grid
