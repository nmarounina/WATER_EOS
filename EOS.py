from parameters import *
import math as m
from scipy.optimize import *
from scipy import constants
import phase_limits

class DataPoint:
    def __init__(self, temperature, pressure):
        self.T = temperature
        self.p = pressure

        self.rho = self.get_density()

        self.grid="None"
        self.phase = "vapor"  # for now


    def get_density(self):

    return 1




    def search_rho_for_given_p(self, rho_search):

    return 1

    def chose_the_right_grid(self):

        datapoint_in_ice_VII=False
        if self.T >= 353.3 and self.T < 800.:
            if phase_limits.get_liq_to_iceVII_phase_line_upto_800K(self.T) < self.p:
                datapoint_in_ice_VII=True

        if self.T >= 800. and self.T < T_IAPWS_to_Mazevet+domega:
            if phase_limits.get_liq_to_iceVII_phase_line_upto_1237K(self.T) < self.p:
                datapoint_in_ice_VII = True




        if (self.T < 500. and self.T > 190. and self.p < pVIVIIL):
            self.grid="seafreeze"
        elif ( self.T>=500 and self.T<= T_IAPWS_to_Mazevet-domega and not datapoint_in_ice_VII ):
            self.grid="IAPWS95"
        elif ( self.T>=500 and self.T<= T_IAPWS_to_Mazevet-domega and datapoint_in_ice_VII ):
            self.grid="HPLT"
        elif ( self.T > T_IAPWS_to_Mazevet-domega and self.T< T_IAPWS_to_Mazevet+domega ) :
            self.grid="transition_btw_IAPWS_and_Mazevet"
        elif self.T<= T_IAPWS_to_Mazevet+domega :
            self.grid="check_rho_first"

