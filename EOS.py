import parameters
from parameters import *
import IAPWS95
import phase_limits
import seafreeze as sf

import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)


class DataPoint:
    def __init__(self, temperature, pressure):
        self.T = temperature
        self.p = pressure

        self.IAPWS95_object = None
        self.seafreeze_object = None
        self.grid = None
        self.choose_the_right_grid()

        self.rho=0
        self.s=0.
        self.u = 0
        self.Cp = 0

        self.get_rho()

        # self.get_s()
        # self.get_u()
        # self.get_Cp()

        self.phase = "vapor"  # for now



    def choose_the_right_grid(self):

        datapoint_in_ice_VII=False
        if self.T >= 353.3 and self.T < 800.:
            if phase_limits.get_liq_to_iceVII_phase_line_upto_800K(self.T) < self.p:
                datapoint_in_ice_VII=True

        if self.T >= 800. and self.T < T_IAPWS_to_Mazevet+domega:
            if phase_limits.get_liq_to_iceVII_phase_line_upto_1237K(self.T) < self.p:
                datapoint_in_ice_VII = True

        datapoint_in_vapor = False
        if self.T >= 190 and self.T < Tt:
            if phase_limits.get_SVP_vap_ice(self.T) > self.p:
                datapoint_in_vapor=True

        if self.T >=Tt and self.T < 500.:
            if phase_limits.get_SVP_vap_liq(self.T) > self.p:
                datapoint_in_vapor = True

        if 500. > self.T > 190. and self.p < pVIVIIL and not datapoint_in_vapor :

            self.grid="seafreeze"
            PT = np.empty((1,), object)
            PT[0] = (self.p / 1e6, self.T)
            whereIam = sf.whichphase( PT )[0]
            material="_IAPWS95"
            if (whereIam == 1):
                material = "Ih"
            elif (whereIam == 2):
                material = "II"
            elif (whereIam == 3):
                material = "III"
            elif (whereIam == 5):
                material = "V"
            elif (whereIam == 6):
                material = "VI"
            elif (whereIam == 0):
                material = "water1"  # _IAPWS95"
            else:
                print("Problem picking an phase for seafreeze module, parameters:",
                      self.p,self.T, whereIam)
                exit()

            self.seafreeze_object = sf.seafreeze( PT , material)

        elif 500. > self.T > 190. and self.p < pVIVIIL and datapoint_in_vapor:
            self.grid = "IAPWS95"
            self.IAPWS95_object = IAPWS95.DataPoint_IAPWS95(self.T, self.p)

        if 500. > self.T > 190. and self.p > pVIVIIL:
            self.grid = "HPLT"

        elif (500 <= self.T <= T_IAPWS_to_Mazevet - domega and not datapoint_in_ice_VII):
            self.grid="IAPWS95"
            self.IAPWS95_object = IAPWS95.DataPoint_IAPWS95(self.T, self.p)

        elif (500 <= self.T <= T_IAPWS_to_Mazevet - domega and datapoint_in_ice_VII):
            self.grid="HPLT"

        elif T_IAPWS_to_Mazevet - domega <= self.T :
            self.IAPWS95_object=IAPWS95.DataPoint_IAPWS95(self.T, self.p)

            if self.IAPWS95_object.rho*M_h2o < 1000.:
                self.grid="IAPWS95"
            else:
                self.grid = "Mazevet"
                if T_IAPWS_to_Mazevet - domega < self.T < T_IAPWS_to_Mazevet + domega and not datapoint_in_ice_VII:
                    self.grid = "transition_btw_IAPWS_and_Mazevet"
                elif T_IAPWS_to_Mazevet - domega < self.T < T_IAPWS_to_Mazevet + domega and datapoint_in_ice_VII:
                    self.grid = "HPLT"
                # elif self.T <= T_IAPWS_to_Mazevet+domega:
                #     self.grid = "Mazevet"


    def get_rho(self):

        if self.grid == "seafreeze" :
            self.rho=self.seafreeze_object.rho/M_h2o #mol.m-3

        elif self.grid == "IAPWS95" :
            self.rho=self.IAPWS95_object.rho

        elif self.grid == "HPLT":
            self.rho=parameters.rho_HPI_pT(self.p, self.T)/ M_h2o

        elif self.grid == "transition_btw_IAPWS_and_Mazevet":
            self.rho=0.

        elif self.grid == "Mazevet":
            self.rho = parameters.rho_mazevet_pT(self.p, self.T)
        return 1

    def get_s(self):
        return 1

    def get_u(self):
        return 1

    def get_Cp(self):
        return 1
