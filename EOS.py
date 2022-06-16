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

        self.phase = "vapor"  # for now


    def get_density(self):

    return 1




    def search_rho_for_given_p(self, rho_search):

    return 1