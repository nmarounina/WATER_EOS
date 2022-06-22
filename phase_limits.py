from parameters import Tc, pc, Tt, transition_width, T_IAPWS_to_Mazevet
import math as m

def get_SVP_vap_liq(T):  #  in K
    #SVP=saturation vapor pressure
    # between the vapor and the liquid phase
    # Wagner and Pruss 2002, valid up to the critical point, with an incertitude <0.1%
    # Eq. 2.5
    th = 1. - T / Tc
    if not (T >= Tt and T <= Tc):
        raise "Outside of the valid temperature range for this saturation vapor pressure"
    else:
        Psat = Tc / T * (
                    -7.85951783 * th +
                    1.84408259 * th ** 1.5 -
                    11.7866497 * th ** 3 +
                    22.6807411 * th ** 3.5 -
                    15.9618719 * th ** 4 +
                    1.80122502 * th ** 7.5)

    return m.exp(Psat) * pc # Pascals


def get_SVP_vap_ice(T):
    # SVP=saturation vapor pressure
    # between the vapor and ice phase
    # valid from 273.16K to 190K, eq.2.21 in Wagner and Pruss 2002

    if (T > Tt):
        raise "Outside of the valid temperature range for this saturation vapor pressure"
    else:
        th = T / Tt
        Psat = (-13.928169 * (1. - th ** (-1.5)) +
                34.7078238 * (1. - th ** (-1.25)) )


    return m.exp(Psat) * 611.657 # Pascals, the pressure multiplying the exp(Psat) is knowingly different from Pt
    # see Wagner and Pruss 2002 and Wagner et al. 1994 for more details


def get_liq_to_iceVI_phase_line(T):

    if not (T >= 273.31 and T < 353.5):  # ice VI
        raise "Outside of the valid temperature range for this melting line (ice VI)"
    else:
        pmelt = (1. - 1.07476 * (1. - (T / 273.31) ** 4.6)) * 632.4  # MPa, Wagner 1994

    return pmelt*1e6 # Pa

def get_liq_to_iceVII_phase_line_upto_800K(T):

    if not (T >= 353.5 and T < 800.):  # ice VII
        raise "Outside of the valid temperature range for this melting line (ice VII, 353.5-800 K)"
    else:
        theta = T / 355.
        pmelt = 0.85 * ((theta) ** 3.47 - 1.) + 2.17  # GPa, Lin+ 2004

    return pmelt*1e9 #Pa


def  get_liq_to_iceVII_phase_line_upto_1237K(T):

    if not(T >= 800 and T <= T_IAPWS_to_Mazevet  + transition_width):
        raise "Outside of the valid temperature range for this melting line (ice VII, 800K-1273K)"
    else:
        pmelt = 3.3488 + 1.8696e-5 * T ** 2  # GPa, fitted from Hernandez+ 2018

    return pmelt*1e9 #Pa

def w(T, p):  # transition function from  IAPWS to Mazevet in the fluid phase

    T0 = 1273.
    contribution_of_IAPWS = 1.
    #on the scale of 0 to 1, 1-contribution_of_IAPWS = contribution of Mazevet to the density values

    if T <= T0 - transition_width:
        contribution_of_IAPWS = 1.
    elif T0 - transition_width < T <= T0 + transition_width:
        a = -1. / (2. * transition_width)
        b = 0.5 + T0 / (2. * transition_width)
        contribution_of_IAPWS = a * T + b
    elif T > T0 + transition_width:
        contribution_of_IAPWS = 0.

    if p > get_liq_to_iceVII_phase_line_upto_1237K(T) :
        contribution_of_IAPWS = 0.

    return contribution_of_IAPWS