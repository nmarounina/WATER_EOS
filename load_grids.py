from parameters import M_h2o
import math as m
import numpy as np
from scipy import interpolate
from scipy import constants as cst


#################################################################
#################################################################
#################################################################
################################################################# PRETRAITMENT FOR THE USE OF MAZEVET EOS:

def create_grid_for_interpolation_on_Mazevet_data():
    Mazevet_raw = np.loadtxt("./dataIN/eosMazevet.dat")  ### this file should be a SQUARE GRID IN T-rho

    # Columns in the Mazevet file:
    # rho[g/cc], T [K], P [GPa], P/(n_i kT), F/(N_i kT), U/(N_i kT), C_V/(N_i k), chi_T, chi_r, S/(N_i k)
    s_offset = -54.280936454849495  # J.mol-1.K-1
    N_lignes_raw = len(Mazevet_raw)

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


    ######################### COPY THOSE FUNCTIONS IN A FUTURE RELEVANT SPOT
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


def load_HPLT_data():
    hpi = np.loadtxt("./dataIN/dataHPices.dat")

    NlgnHPI = len(hpi)

    THPI = np.zeros(NlgnHPI)
    PHPI = np.zeros(NlgnHPI)
    rhoHPI = np.zeros(NlgnHPI)
    alphaHPI = np.zeros(NlgnHPI)
    CpHPI = np.zeros(NlgnHPI)
    sHPI = np.zeros(NlgnHPI)

    aa_HPI = np.zeros((int(NlgnHPI ** 0.5), int(NlgnHPI ** 0.5)))
    ss_HPI = np.zeros((int(NlgnHPI ** 0.5), int(NlgnHPI ** 0.5)))
    cc_HPI = np.zeros((int(NlgnHPI ** 0.5), int(NlgnHPI ** 0.5)))
    rrho_HPI = np.zeros((int(NlgnHPI ** 0.5), int(NlgnHPI ** 0.5)))

    tt = np.zeros((int(NlgnHPI ** 0.5), int(NlgnHPI ** 0.5)))
    pp = np.zeros((int(NlgnHPI ** 0.5), int(NlgnHPI ** 0.5)))

    for i in range(0, NlgnHPI):
        THPI[i] = hpi[i, 0]
        PHPI[i] = hpi[i, 1]
        rhoHPI[i] = hpi[i, 2]
        alphaHPI[i] = hpi[i, 3]
        CpHPI[i] = hpi[i, 4]
        sHPI[i] = hpi[i, 5]

    thpi = np.linspace(np.min(THPI), np.max(THPI), int(NlgnHPI ** 0.5))
    phpi = np.linspace(m.log10(np.min(PHPI)), m.log10(np.max(PHPI)), int(NlgnHPI ** 0.5))
    phpi = 10 ** phpi

    tt, pp = np.meshgrid(thpi, phpi)  # to double-check the indices

    cpt = 0
    for i in range(0, int(NlgnHPI ** 0.5)):
        for j in range(0, int(NlgnHPI ** 0.5)):
            aa_HPI[i, j] = alphaHPI[cpt]
            ss_HPI[i, j] = sHPI[cpt]
            cc_HPI[i, j] = CpHPI[cpt]
            rrho_HPI[i, j] = rhoHPI[cpt]

            # print(i,j,cpt, )

            cpt += 1

    s_HPI_Tp = interpolate.RectBivariateSpline(phpi, thpi, ss_HPI)
    alpha_HPI_Tp = interpolate.RectBivariateSpline(phpi, thpi, aa_HPI)
    rho_HPI_Tp = interpolate.RectBivariateSpline(phpi, thpi, rrho_HPI)
    Cp_HPI_Tp = interpolate.RectBivariateSpline(phpi, thpi, cc_HPI)

    return 1


#####################################

def load_some_other_data():  # please rename and wtf are those grids
    datarholiq = np.loadtxt("./dataIN/liqgrid_500.dat")
    datarhovap = np.loadtxt("./dataIN/vapgrid_500.dat")
    datarhosc = np.loadtxt("./dataIN/scgrid_500.dat")

    tliq = np.zeros(len(datarholiq))
    plgrid = np.zeros(len(datarholiq))
    rholiq = np.zeros(len(datarholiq))

    tvap = np.zeros(len(datarhovap))
    pvgrid = np.zeros(len(datarhovap))
    rhovap = np.zeros(len(datarhovap))

    tsc = np.zeros(len(datarhosc))
    psc = np.zeros(len(datarhosc))
    rhosc = np.zeros(len(datarhosc))

    for i in range(0, len(datarhovap)):
        tvap[i] = m.log10(datarhovap[i, 0])
        pvgrid[i] = datarhovap[i, 1]
        rhovap[i] = m.log10(datarhovap[i, 2])

    Grmax = np.max(pvgrid)
    Grmin = np.min(pvgrid)
    Tmax = 10 ** np.max(tvap)
    Tmin = 10 ** np.min(tvap)
    tvaptest = np.linspace(m.log10(Tmin), m.log10(Tmax), int(len(datarhovap) ** 0.5))
    lgridtest = np.linspace(Grmin, Grmax, int(len(datarhovap) ** 0.5))
    rhovaptest = np.zeros((int(len(datarhovap) ** 0.5), int(len(datarhovap) ** 0.5)))

    cpt = 0
    for i in range(0, int(len(datarhovap) ** 0.5)):
        for j in range(0, int(len(datarhovap) ** 0.5)):
            rhovaptest[i, j] = m.log10(datarhovap[cpt, 2])
            cpt += 1

    testvap_RBS = interpolate.RectBivariateSpline(tvaptest, lgridtest, rhovaptest)

    ##--------------------------------------------------------------------------------------------
    ##---------------Now the liquid:
    ##--------------------------------------------------------------------------------------------
    ##--------------------------------------------------------------------------------------------
    ##--------------------------------------------------------------------------------------------

    for i in range(0, len(datarholiq)):
        tliq[i] = m.log10(datarholiq[i, 0])
        plgrid[i] = datarholiq[i, 1]
        rholiq[i] = m.log10(datarholiq[i, 2])

    Grmax = np.max(plgrid)
    Grmin = np.min(plgrid)
    Tmax = 10 ** np.max(tliq)
    Tmin = 10 ** np.min(tliq)

    tliqtest = np.linspace(m.log10(Tmin), m.log10(Tmax), int(len(datarholiq) ** 0.5))
    plgridtest = np.linspace(Grmin, Grmax, int(len(datarholiq) ** 0.5))
    rholiqtest = np.zeros((int(len(datarholiq) ** 0.5), int(len(datarholiq) ** 0.5)))

    cpt = 0
    for i in range(0, int(len(datarholiq) ** 0.5)):
        for j in range(0, int(len(datarholiq) ** 0.5)):
            rholiqtest[i, j] = m.log10(datarholiq[cpt, 2])
            cpt += 1

    testliq_RBS = interpolate.RectBivariateSpline(tliqtest, plgridtest, rholiqtest)

    ##--------------------------------------------------------------------------------------------
    ##---------------And finally the supercritical grid:
    ##--------------------------------------------------------------------------------------------
    ##--------------------------------------------------------------------------------------------
    ##--------------------------------------------------------------------------------------------

    for i in range(0, len(datarhosc)):
        tsc[i] = m.log10(datarhosc[i, 0])
        psc[i] = m.log10(datarhosc[i, 1])
        rhosc[i] = m.log10(datarhosc[i, 2])

    Pmax = 10 ** np.max(psc)
    Pmin = 10 ** np.min(psc)
    Tmax = 10 ** np.max(tsc)
    Tmin = 10 ** np.min(tsc)

    tsctest = np.linspace(m.log10(Tmin), m.log10(Tmax), int(len(datarhosc) ** 0.5))
    psctest = np.linspace(m.log10(Pmin), m.log10(Pmax), int(len(datarhosc) ** 0.5))
    rhosctest = np.zeros((int(len(datarhosc) ** 0.5), int(len(datarhosc) ** 0.5)))

    cpt = 0
    for i in range(0, int(len(datarhosc) ** 0.5)):
        for j in range(0, int(len(datarhosc) ** 0.5)):
            rhosctest[i, j] = m.log10(datarhosc[cpt, 2])
            cpt += 1

    testsc_RBS = interpolate.RectBivariateSpline(tsctest, psctest, rhosctest)

    return 1
