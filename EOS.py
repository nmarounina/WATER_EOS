from parameters import *
import IAPWS95
import phase_limits

class DataPoint:
    def __init__(self, temperature, pressure):
        self.T = temperature
        self.p = pressure

        self.grid = self.chose_the_right_grid()

        self.rho=self.get_rho()
        self.s=self.get_s()
        self.u = self.get_u()
        self.Cp = self.get_Cp()

        self.phase = "vapor"  # for now



    def chose_the_right_grid(self):

        datapoint_in_ice_VII=False
        if self.T >= 353.3 and self.T < 800.:
            if phase_limits.get_liq_to_iceVII_phase_line_upto_800K(self.T) < self.p:
                datapoint_in_ice_VII=True

        if self.T >= 800. and self.T < T_IAPWS_to_Mazevet+domega:
            if phase_limits.get_liq_to_iceVII_phase_line_upto_1237K(self.T) < self.p:
                datapoint_in_ice_VII = True




        if (500. > self.T > 190. and self.p < pVIVIIL):
            self.grid="seafreeze"
        elif (500 <= self.T <= T_IAPWS_to_Mazevet - domega and not datapoint_in_ice_VII):
            self.grid="IAPWS95"
        elif (500 <= self.T <= T_IAPWS_to_Mazevet - domega and datapoint_in_ice_VII):
            self.grid="HPLT"

        elif T_IAPWS_to_Mazevet - domega < self.T :
            Point_from_IAPWS=IAPWS95.DataPoint_IAPWS95(self.T, self.p)

            if Point_from_IAPWS.rho/M_h2o < 1000.:
                self.grid="IAPWS95"
            else:

                if T_IAPWS_to_Mazevet - domega < self.T < T_IAPWS_to_Mazevet + domega and not datapoint_in_ice_VII:
                    self.grid = "transition_btw_IAPWS_and_Mazevet"
                elif T_IAPWS_to_Mazevet - domega < self.T < T_IAPWS_to_Mazevet + domega and datapoint_in_ice_VII:
                    self.grid = "transition_btw_HPLT_and_Mazevet"
                elif self.T <= T_IAPWS_to_Mazevet+domega:
                    self.grid = "Mazevet"

    def get_rho(self):
        return 1

    def get_s(self):
        return 1

    def get_u(self):
        return 1

    def get_Cp(self):
        return 1
