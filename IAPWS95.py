from parameters import *
import math as m
from scipy.optimize import *
from scipy import constants
import phase_limits

class DataPoint_IAPWS95:
    def __init__(self, temperature, pressure):
        self.T = temperature
        self.p = pressure

        # Parameters of the EOS
        self.delt = 0.
        self.tau = Tc / self.T

        self.rho = self.get_density()

        self.delt = self.rho / rhoc

        # list of derivatives
        self.phi0 = self.compute_phi0()
        self.phir = self.compute_phir()
        self.dphirt = self.compute_dphir_t()
        self.dphirdd = self.compute_dphir_dd()
        self.dphirdt = self.compute_dphir_dt()
        self.dphird = self.compute_dphir_d()
        self.dphi0t = self.compute_dphi0_t()
        self.dphi0tt = self.compute_dphi0_tt()
        self.dphirtt = self.compute_dphir_tt()

        self.Cp = (-1. * self.tau ** 2. * (self.dphi0tt + self.dphirtt) +
                   (1. + self.delt * self.dphird - self.delt * self.tau * self.dphirdt) ** 2. /
                   (1. + 2. * self.delt * self.dphird + self.delt ** 2. * self.dphirdd)) * constants.R

        self.s = (self.tau * (self.dphi0t + self.dphirt) - self.phi0 - self.phir) * constants.R
        self.u = (self.tau * (self.dphi0t + self.dphirt)) * constants.R * self.T
        self.dhdp_T = (1. - (1. + self.delt * self.dphird - self.delt * self.tau * self.dphirdt) / (
                1. + 2. * self.delt * self.dphird + self.delt ** 2 * self.dphirdd)) * 1. / self.rho
        self.dsdp = (1. / (self.T * self.rho) - 1. / self.T * self.dhdp_T) * (-1.)

        self.phase = "vapor"  # for now

    def get_density(self):

        if self.T > Tc :

            return bisect(self.search_rho_for_given_p,
                   1e-10,
                   self.p / (self.T * constants.R) *10.,
                   )  # density, mol.m-3

        elif self.T<=Tc and self.T>=Tt :

            SVP=phase_limits.get_SVP_vap_liq(self.T)

            if self.p>SVP:
                return newton(self.search_rho_for_given_p,
                              2000./M_h2o,
                              )  # density, mol.m-3
            else:
                return newton(self.search_rho_for_given_p,
                              self.p / (self.T * constants.R),
                              )  # density, mol.m-3

        elif self.T < Tt :

            SVP = phase_limits.get_SVP_vap_ice(self.T)

            if self.p > SVP:
                return newton(self.search_rho_for_given_p,
                              2000. / M_h2o,
                              )  # density, mol.m-3
            else:
                return newton(self.search_rho_for_given_p,
                              self.p / (self.T * constants.R),
                              )  # density, mol.m-3




    def search_rho_for_given_p(self, rho_search):
        self.delt = rho_search / rhoc
        return ((1. + self.delt * self.compute_dphir_d()) *
                rho_search * constants.R * self.T - self.p) / self.p

    #################################################################
    #################################################################
    #################################################################
    #################################################################

    # def pressure(self, rhoi, T):
    #     delt = rhoi / rhoc
    #     tau = Tc / T
    #
    #     phir, dphird, dphirt, dphirtt, dphirdt, dphirdd = compute_phir(delt, tau)
    #
    #     return (1. + delt * dphird) * rhoi * Rig * T

    #################################################################
    #################################################################
    #################################################################
    #################################################################

    def compute_phi0(self):

        sum = 0.

        for i in range(3, 8):
            sum = sum + n0[i] * m.log(1. - m.exp(-1. * gam0[i] * self.tau))

        phi0 = m.log(self.delt) + n0[0] + n0[1] * self.tau + n0[2] * m.log(self.tau) + sum

        return phi0

    def compute_dphi0_t(self):
        # derivative of the ideal part with tau
        dsum = 0.
        for i in range(3, 8):
            dsum = dsum + n0[i] * gam0[i] * ((1. - m.exp(-1. * gam0[i] * self.tau)) ** (-1.) - 1.)

        return n0[1] + n0[2] / self.tau + dsum

    def compute_dphi0_tt(self):
        # second derivative of the ideal part of the eos, with tau
        dstt = 0.
        for i in range(3, 8):
            dstt = dstt + n0[i] * gam0[i] ** 2. * m.exp(-1. * gam0[i] * self.tau) * (
                    1. - m.exp(-1. * gam0[i] * self.tau)) ** (-2.)

        return -1. * n0[2] / self.tau ** 2. - dstt

    def compute_phir(self):

        sum1 = 0.
        for i in range(0, 7):
            sum1 = sum1 + n[i] * self.delt ** d[i] * self.tau ** ti[i]

        sum2 = 0.
        for i in range(7, 51):
            sum2 = sum2 + n[i] * self.delt ** (d[i]) * self.tau ** ti[i] * m.exp(-1. * self.delt ** c[i])

        sum3 = 0.
        for i in range(51, 54):
            sum3 = sum3 + n[i] * self.delt ** d[i] * self.tau ** ti[i] * m.exp(
                -1. * al[i] * (self.delt - eps[i]) ** 2. - be[i] * (self.tau - gam[i]) ** 2.)

        sum4 = 0.
        for i in range(54, 56):
            psi = m.exp(-1. * C[i] * (self.delt - 1.) ** 2. - D[i] * (self.tau - 1.) ** 2.)
            theta = (1. - self.tau) + A[i] * ((self.delt - 1.) ** 2.) ** (1. / (2. * bet[i]))
            Delta = theta ** 2 + B[i] * ((self.delt - 1.) ** 2.) ** (a[i])

            sum4 = sum4 + n[i] * Delta ** (b[i]) * self.delt * psi

        return sum1 + sum2 + sum3 + sum4

    def compute_dphir_d(self):

        ds1d = 0.
        for i in range(0, 7):
            ds1d = ds1d + n[i] * d[i] * self.delt ** (d[i] - 1.) * self.tau ** ti[i]

        ds2d = 0.
        for i in range(7, 51):
            ds2d = ds2d + n[i] * m.exp(-1. * self.delt ** c[i]) * (
                    self.delt ** (d[i] - 1.) * self.tau ** ti[i] * (d[i] - c[i] * self.delt ** c[i]))

        ds3d = 0.
        for i in range(51, 54):
            ds3d = ds3d + n[i] * self.delt ** d[i] * self.tau ** ti[i] * m.exp(
                -1. * al[i] * (self.delt - eps[i]) ** 2. - be[i] * (self.tau - gam[i]) ** 2.) * (
                           d[i] / self.delt - 2. * al[i] * (self.delt - eps[i]))

        ds4d = 0.
        for i in range(54, 56):
            psi = m.exp(-1. * C[i] * (self.delt - 1.) ** 2. - D[i] * (self.tau - 1.) ** 2.)
            theta = (1. - self.tau) + A[i] * ((self.delt - 1.) ** 2.) ** (1. / (2. * bet[i]))
            Delta = theta ** 2 + B[i] * ((self.delt - 1.) ** 2.) ** (a[i])
            dpsid = -2. * C[i] * (self.delt - 1.) * psi

            dDeltad = (self.delt - 1.) * (
                    A[i] * theta * 2. / bet[i] * ((self.delt - 1.) ** 2.) ** (1. / (2. * bet[i]) - 1.) + 2. * B[i] * a[
                i] * (
                            (self.delt - 1.) ** 2.) ** (a[i] - 1.))

            dDeltbid = b[i] * Delta ** (b[i] - 1.) * dDeltad
            ds4d = ds4d + n[i] * (Delta ** b[i] * (psi + self.delt * dpsid) + dDeltbid * self.delt * psi)

        return ds1d + ds2d + ds3d + ds4d

    def compute_dphir_t(self):
        ds1t = 0.
        for i in range(0, 7):
            ds1t = ds1t + n[i] * ti[i] * self.delt ** d[i] * self.tau ** (ti[i] - 1.)

        ds2t = 0.
        for i in range(7, 51):
            ds2t = ds2t + n[i] * ti[i] * self.delt ** d[i] * self.tau ** (ti[i] - 1.) * m.exp(-1. * self.delt ** c[i])

        ds3t = 0.
        for i in range(51, 54):
            ds3t = ds3t + n[i] * self.delt ** d[i] * self.tau ** ti[i] * m.exp(
                -1. * al[i] * (self.delt - eps[i]) ** 2. - be[i] * (self.tau - gam[i]) ** 2.) * (
                           ti[i] / self.tau - 2. * be[i] * (self.tau - gam[i]))

        ds4t = 0.
        for i in range(54, 56):
            psi = m.exp(-1. * C[i] * (self.delt - 1.) ** 2. - D[i] * (self.tau - 1.) ** 2.)
            theta = (1. - self.tau) + A[i] * ((self.delt - 1.) ** 2.) ** (1. / (2. * bet[i]))
            Delta = theta ** 2 + B[i] * ((self.delt - 1.) ** 2.) ** (a[i])
            dDeltbit = -2. * theta * b[i] * Delta ** (b[i] - 1.)
            dpsit = -2. * D[i] * (self.tau - 1.) * psi

            ds4t = ds4t + n[i] * self.delt * (dDeltbit * psi + Delta ** b[i] * dpsit)

        return ds1t + ds2t + ds3t + ds4t

    def compute_dphir_tt(self):
        ds1tt = 0.

        for i in range(0, 7):
            ds1tt = ds1tt + n[i] * ti[i] * (ti[i] - 1.) * self.delt ** d[i] * self.tau ** (ti[i] - 2.)

        ds2tt = 0.
        for i in range(7, 51):
            ds2tt = ds2tt + n[i] * ti[i] * (ti[i] - 1.) * self.delt ** d[i] * self.tau ** (ti[i] - 2.) * m.exp(
                -1. * self.delt ** c[i])

        ds3tt = 0.
        for i in range(51, 54):
            ds3tt = ds3tt + n[i] * self.delt ** d[i] * self.tau ** ti[i] * m.exp(
                -1. * al[i] * (self.delt - eps[i]) ** 2. - be[i] * (self.tau - gam[i]) ** 2.) * (
                            (ti[i] / self.tau - 2. * be[i] * (self.tau - gam[i])) ** 2. - ti[i] / self.tau ** 2. - 2. *
                            be[i])

        ds4tt = 0.
        for i in range(54, 56):
            psi = m.exp(-1. * C[i] * (self.delt - 1.) ** 2. - D[i] * (self.tau - 1.) ** 2.)
            theta = (1. - self.tau) + A[i] * ((self.delt - 1.) ** 2.) ** (1. / (2. * bet[i]))
            Delta = theta ** 2 + B[i] * ((self.delt - 1.) ** 2.) ** (a[i])

            dDeltbit = -2. * theta * b[i] * Delta ** (b[i] - 1.)
            dpsit = -2. * D[i] * (self.tau - 1.) * psi
            dpsitt = (2. * D[i] * (self.tau - 1.) ** 2. - 1.) * 2. * D[i] * psi
            dDelbitt = 2. * b[i] * Delta ** (b[i] - 1.) + 4. * theta ** 2 * b[i] * (b[i] - 1.) * Delta ** (b[i] - 2.)

            ds4tt = ds4tt + n[i] * self.delt * (dDelbitt * psi + 2. * dDeltbit * dpsit + Delta ** b[i] * dpsitt)

        return ds1tt + ds2tt + ds3tt + ds4tt

    def compute_dphir_dd(self):

        ds1dd = 0.

        for i in range(0, 7):
            ds1dd = ds1dd + n[i] * d[i] * (d[i] - 1.) * self.delt ** (d[i] - 2.) * self.tau ** ti[i]

        ds2dd = 0.
        for i in range(7, 51):
            ds2dd = ds2dd + n[i] * m.exp(-1. * self.delt ** c[i]) * (self.delt ** (d[i] - 2.) * self.tau ** ti[i] * (
                    (d[i] - c[i] * self.delt ** c[i]) * (d[i] - 1. - c[i] * self.delt ** c[i]) - c[
                i] ** 2. * self.delt ** c[i]))

        ds3dd = 0.
        for i in range(51, 54):
            ds3dd = ds3dd + n[i] * self.tau ** ti[i] * m.exp(
                -1. * al[i] * (self.delt - eps[i]) ** 2. - be[i] * (self.tau - gam[i]) ** 2.) * (
                            -2. * al[i] * self.delt ** d[i] + 4. * al[i] ** 2 * self.delt ** d[i] * (
                            self.delt - eps[i]) ** 2. - 4. *
                            d[i] * al[i] * self.delt ** (d[i] - 1.) * (self.delt - eps[i]) + d[i] * (
                                    d[i] - 1.) * self.delt ** (
                                    d[i] - 2.))

        ds4dd = 0.
        for i in range(54, 56):
            psi = m.exp(-1. * C[i] * (self.delt - 1.) ** 2. - D[i] * (self.tau - 1.) ** 2.)
            theta = (1. - self.tau) + A[i] * ((self.delt - 1.) ** 2.) ** (1. / (2. * bet[i]))
            Delta = theta ** 2 + B[i] * ((self.delt - 1.) ** 2.) ** (a[i])

            dDeltad = (self.delt - 1.) * (
                    A[i] * theta * 2. / bet[i] * ((self.delt - 1.) ** 2.) ** (1. / (2. * bet[i]) - 1.) + 2. * B[i] * a[
                i] * (
                            (self.delt - 1.) ** 2.) ** (a[i] - 1.))
            dDeltbid = b[i] * Delta ** (b[i] - 1.) * dDeltad
            dpsid = -2. * C[i] * (self.delt - 1.) * psi
            dpsidd = (2. * C[i] * (self.delt - 1) ** 2. - 1.) * 2. * C[i] * psi
            dDeldd = 1. / (self.delt - 1.) * dDeltad + (self.delt - 1.) ** 2 * (
                    4. * B[i] * a[i] * (a[i] - 1.) * ((self.delt - 1.) ** 2.) ** (a[i] - 2.) + 2. * A[i] ** 2. * (
                    1. / bet[i]) ** 2. * (((self.delt - 1.) ** 2.) ** (1. / (2. * bet[i]) - 1.)) ** 2 + A[
                        i] * theta * 4. / bet[i] * (1. / (2. * bet[i]) - 1.) * ((self.delt - 1.) ** 2.) ** (
                            1. / (2. * bet[i]) - 2.))
            dDelbidd = b[i] * (Delta ** (b[i] - 1.) * dDeldd + (b[i] - 1.) * Delta ** (b[i] - 2.) * dDeltad ** 2.)

            ds4dd = ds4dd + n[i] * (Delta ** b[i] * (2. * dpsid + self.delt * dpsidd) + 2. * dDeltbid * (
                    psi + self.delt * dpsid) + dDelbidd * self.delt * psi)

        return ds1dd + ds2dd + ds3dd + ds4dd

    def compute_dphir_dt(self):
        ds1dt = 0.

        for i in range(0, 7):
            ds1dt = ds1dt + n[i] * d[i] * ti[i] * self.tau ** (ti[i] - 1.) * self.delt ** (d[i] - 1.)

        ds2dt = 0.

        for i in range(7, 51):
            ds2dt = ds2dt + n[i] * ti[i] * self.tau ** (ti[i] - 1.) * self.delt ** (d[i] - 1.) * (
                    d[i] - c[i] * self.delt ** c[i]) * m.exp(
                -1. * self.delt ** c[i])

        ds3dt = 0.
        for i in range(51, 54):
            ds3dt = ds3dt + n[i] * self.delt ** d[i] * self.tau ** ti[i] * m.exp(
                -1. * al[i] * (self.delt - eps[i]) ** 2. - be[i] * (self.tau - gam[i]) ** 2.) * (
                            ti[i] / self.tau - 2. * be[i] * (self.tau - gam[i])) * (
                            d[i] / self.delt - 2. * al[i] * (self.delt - eps[i]))

        ds4dt = 0.
        for i in range(54, 56):
            psi = m.exp(-1. * C[i] * (self.delt - 1.) ** 2. - D[i] * (self.tau - 1.) ** 2.)
            theta = (1. - self.tau) + A[i] * ((self.delt - 1.) ** 2.) ** (1. / (2. * bet[i]))
            Delta = theta ** 2 + B[i] * ((self.delt - 1.) ** 2.) ** (a[i])

            dDeltad = (self.delt - 1.) * (
                    A[i] * theta * 2. / bet[i] * ((self.delt - 1.) ** 2.) ** (1. / (2. * bet[i]) - 1.) + 2. * B[i] * a[
                i] * (
                            (self.delt - 1.) ** 2.) ** (a[i] - 1.))
            dDeltbid = b[i] * Delta ** (b[i] - 1.) * dDeltad
            dDeltbit = -2. * theta * b[i] * Delta ** (b[i] - 1.)
            dpsid = -2. * C[i] * (self.delt - 1.) * psi
            dpsit = -2. * D[i] * (self.tau - 1.) * psi

            dpsidt = 4. * C[i] * D[i] * (self.delt - 1.) * (self.tau - 1.) * psi

            dDelbidt = -1. * A[i] * b[i] * 2. / bet[i] * Delta ** (b[i] - 1.) * (self.delt - 1.) * (
                    (self.delt - 1.) ** 2.) ** (
                               1. / (2. * bet[i]) - 1.) - 2. * theta * b[i] * (b[i] - 1.) * Delta ** (
                               b[i] - 2.) * dDeltad

            ds4dt = ds4dt + n[i] * (
                    Delta ** b[i] * (dpsit + self.delt * dpsidt) + self.delt * dDeltbid * dpsit + dDeltbit * (
                    psi + self.delt * dpsid) + dDelbidt * self.delt * psi)

        return ds1dt + ds2dt + ds3dt + ds4dt
