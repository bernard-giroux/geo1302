#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 19:18:41 2017

@author: giroux
"""
import time
import numpy as np
from scipy.interpolate import interpn
import matplotlib.pyplot as plt
from numba import jit

from utils import nargout


class GrilleFDTD:
    def __init__(self, x, z):
        """
        Input
            x: coordonnées des noeuds selon x     (mètres)
            z: coordonnées des noeuds selon z     (mètres)
        """
        self.x = x
        self.z = z

    @property
    def x(self):
        "Coordonnées des noeuds selon x"
        return self._x

    @x.setter
    def x(self, val):
        try:
            tmp = np.array(val, dtype=np.float64)
            if tmp.ndim != 1:
                raise ValueError('1D array needed')
            if len(np.unique(np.diff(tmp))) > 1:
                raise ValueError('Constant step size needed')
        except Exception as e:
            raise e

        self._x = tmp
        self.dx = tmp[1]-tmp[0]
        self.nx = tmp.size
        if '_z' in self.__dict__:
            dz = self._z[1]-self._z[0]
            if dz != self.dx:
                raise ValueError('dx should be equal to dz')

    @property
    def z(self):
        "Coordonnées des noeuds selon z"
        return self._z

    @z.setter
    def z(self, val):
        try:
            tmp = np.array(val, dtype=np.float64)
            if tmp.ndim != 1:
                raise ValueError('1D array needed')
            if len(np.unique(np.diff(tmp))) > 1:
                raise ValueError('Constant step size needed')
        except Exception as e:
            raise e

        self._z = tmp
        self.nz = tmp.size
        if 'dx' in self.__dict__:
            dz = self._z[1]-self._z[0]
            if dz != self.dx:
                raise ValueError('dx should be equal to dz')

    def defProp(self, Vp, Vs, rho, dt):
        """
        Calcul des coefficients constants

        Input
            Vp: vitesse des ondes P en m/s                (taille: nz par nx)
            Vs: vitesse des ondes S en m/s                (taille: nz par nx)
            rho: densité en kg/m^3                        (taille: nz par nx)
            dt: pas temporel en s
        """
        if Vp.shape != (self.nz, self.nx):
            raise ValueError('Size of Vp incorrect')
        if Vs.shape != (self.nz, self.nx):
            raise ValueError('Size of Vs incorrect')
        if rho.shape != (self.nz, self.nx):
            raise ValueError('Size of rho incorrect')

        self.dt = dt
        self.b1 =

        # moyenne arithmétique pour la densité
        # à i+1/2, j+1/2

        x =
        z =
        x = x.reshape(-1, 1)
        z = z.reshape(-1, 1)
        # cannot extrapolate with interpn, we duplicate last column
        x[-1] = x[-2]
        z[-1] = z[-2]     # idem for last row
        # create array of node coordinates at i+1/2, j+1/2
        xi = np.hstack((

        rho2 = interpn((self.z, self.x), rho, xi)
        rho2 = rho2.reshape(self.nz, self.nx)

        self.b2 =

        # moyenne harmonique pour les constantes d'élasticité

        # à i+1/2, j
        mu = rho * Vs*Vs
        lambd = rho * Vp*Vp - 2.0*mu
        z = self.z.reshape(-1, 1)

        xi = np.hstack((

        lambda2 =
        mu2 =
        lambda2 = lambda2.reshape(self.nz, self.nx)
        mu2 = mu2.reshape(self.nz, self.nx)

        self.lm = dt / self.dx * (lambda2 + 2*mu2)
        self.l = dt / self.dx * lambda2

        # à i, j+1/2
        x = self.x
        z = self.z + self.dx/2.0
        x = x.reshape(-1, 1)
        z = z.reshape(-1, 1)
        z[-1] = z[-2]

        xi = np.hstack((

        mu2 =
        mu2 = mu2.reshape(self.nz, self.nx)

        self.m = dt / self.dx * mu2

    @jit
    def propage0(self, src, t, showPlot=False):

        nstep = 1 + int(t/self.dt)

        tau_xx = np.zeros((self.nz, self.nx))
        tau_zz = np.zeros((self.nz, self.nx))
        tau_xz = np.zeros((self.nz, self.nx))
        v_x = np.zeros((self.nz, self.nx))
        v_z = np.zeros((self.nz, self.nx))

        if showPlot:
            fig, (ax1, ax2) = plt.subplots(ncols=2)
            im1 = ax1.imshow(v_x)
            im2 = ax2.imshow(v_z)
            stitle = fig.suptitle('')
            plt.show(block=False)

        for m in range(nstep):

            # applique la source au noeud (i+1/2,j)
            tau_zz[src.j, src.i] += src(m)
            tau_xx[src.j, src.i] += src(m)

            for j in np.arange(1, self.nz):
                for i in np.arange(1, self.nx):
                    v_x[j,i] = v_x[j,i] + self.b1[j,i] * (
                        (tau_xx[j,i]-tau_xx[j,i-1]) +
                        (tau_xz[j,i]-tau_xz[j-1,i]))

            for j in np.arange(0, self.nz-1):
                for i in np.arange(0, self.nx-1):
                    v_z[j,i] = v_z[j,i] + self.b2[j,i] * (
                        (tau_xz[j,i+1]-tau_xz[j,i]) +
                        (tau_zz[j+1,i]-tau_zz[j,i]))

            for j in np.arange(1, self.nz):
                for i in np.arange(0, self.nx-1):
                    tau_xx[j,i] += self.lm[j,i] * (v_x[j,i+1]-v_x[j,i]) + \
                        self.l[j,i] * (v_z[j,i]-v_z[j-1,i])

                    tau_zz[j,i] += self.lm[j,i] * (v_z[j,i]-v_z[j-1,i]) + \
                        self.l[j,i] * (v_x[j,i+1]-v_x[j,i])

            for j in np.arange(0, self.nz-1):
                for i in np.arange(1, self.nx):
                    tau_xz[j,i] += self.m[j,i] * (
                        v_x[j+1,i]-v_x[j,i] +
                        v_z[j,i]-v_z[j,i-1])

            if showPlot and np.remainder(m, 20) == 0:
                im1.set_data(v_x)
                im1.set_clim(v_x.min(), v_x.max())
                im2.set_data(v_z)
                im2.set_clim(v_z.min(), v_z.max())
                fig.canvas.draw()
                stitle.set_text('t = {0:6.3f} (it no {1:d}/{2:d})'.format(m*self.dt, m+1, nstep))
                plt.pause(0.01)

    def propage(self, src, t, showPlot=False):

        nstep = 1 + int(t/self.dt)

        tau_xx = np.zeros((self.nz, self.nx))
        tau_zz = np.zeros((self.nz, self.nx))
        tau_xz = np.zeros((self.nz, self.nx))
        v_x = np.zeros((self.nz, self.nx))
        v_z = np.zeros((self.nz, self.nx))

        if showPlot:
            fig, (ax1, ax2) = plt.subplots(ncols=2)
            im1 = ax1.imshow(v_x)
            im2 = ax2.imshow(v_z)
            stitle = fig.suptitle('')
            plt.show(block=False)

        for m in range(nstep):

            # applique la source au noeud (i+1/2,j)
            tau_zz[src.j, src.i] += src(m)
            tau_xx[src.j, src.i] += src(m)

            v_x[] +=

            v_z[] +=

            tau_xx[] +=

            tau_zz[] +=

            tau_xz[] +=

            if showPlot and np.remainder(m, 20) == 0:
                im1.set_data(v_x)
                im1.set_clim(v_x.min(), v_x.max())
                im2.set_data(v_z)
                im2.set_clim(v_z.min(), v_z.max())
                fig.canvas.draw()
                stitle.set_text('t = {0:6.3f} (it no {1:d}/{2:d})'.\
                                format(m*self.dt, m+1, nstep))
                plt.pause(0.01)


class Source:
    """
    Classe mère pour les fonctions sources
    """
    def __init__(self, i, j, A):
        """
        Constructeur

        Paramètres
            i : indice du noeud selon x ou est appliquée la source
            j : indice du noeud selon z ou est appliquée la source
            A : amplitude de la source
        """
        self.i = i
        self.j = j
        self.A = A

    def __call__(self, ind):
        """
        Retourne la valeur de la fonction source à l'indice temporel ind
        """
        if ind < self.f.size:
            return self.A * self.f[ind]
        else:
            return 0.0

    def plot(self):
        fig, ax = plt.subplots(ncols=2)
        t = self.t - self.t[0]  # make sure t starts at 0 (e.g. Ricker)
        ax[0].plot(t, self.A * self.f)
        ax[0].set_xlabel('Temps (s)')

        dt = self.t[1]-self.t[0]
        f = np.hstack((self.f, np.zeros((1024,))))   # zero padding
        S = np.abs(np.fft.fft(f))**2
        f = np.fft.fftfreq(S.size, dt)
        idx = f >= 0                                  # freq positives
        ax[1].plot(f[idx],S[idx])
        ax[1].set_xlabel('Fréquence (Hz)')
        plt.show()


class Impulsion(Source):
    def __init__(self, i, j, alpha, t0, dt, A=1.0):
        Source.__init__(self, i, j, A)
        self.t = np.arange(0.0, 2.0*t0 + dt/3, dt)
        self.f = np.exp(-alpha * (self.t-t0)**2)


class DeriveeImpulsion(Source):
    def __init__(self, i, j, alpha, t0, dt, A=1.0):
        Source.__init__(self, i, j, A)
        self.t = np.arange(0.0, 2.0*t0 + dt/3, dt)
        self.f = -2.0 * alpha * (self.t-t0) * np.exp(-alpha * (self.t-t0)**2)


class Ricker(Source):
    def __init__(self, i, j, fdom, dt, A=1.0):
        Source.__init__(self, i, j, A)
        self.t = np.arange(-1.0/fdom, 1.0/fdom + dt/3, dt)
        self.f =


if __name__ == '__main__':

    checkProp = False
    timeForLoop = True

    if checkProp:
        dx = 50.0
        x = np.arange(0.0, 200.1, dx)
        z = np.arange(0.0, 150.1, dx)

        g = GrilleFDTD(x, z)

        Vp = 4000.0 + np.zeros((z.size, x.size))
        Vp[1, 1] = 5000.0
        Vp[1, 2] = 3000.0
        sigma = 0.25              # coeff Poisson
        Vs = Vp * np.sqrt((0.5-sigma)/(1.0-sigma))
        rho = 2670.0 + np.zeros(Vp.shape)
        rho[1, 1] = 2500.0
        rho[2, 2] = 2700.0
        dt = 0.006

        g.defProp(Vp, Vs, rho, dt)

        plt.figure()
        plt.imshow(g.b1)
        plt.title('b1')
        plt.gca().set_aspect('equal')
        plt.gca().autoscale(tight=True)
        plt.colorbar()
        plt.show()

        plt.figure()
        plt.imshow(g.b2)
        plt.title('b2')
        plt.gca().set_aspect('equal')
        plt.gca().autoscale(tight=True)
        plt.colorbar()
        plt.show()

        plt.figure()
        plt.imshow(g.lm)
        plt.title('lm')
        plt.gca().set_aspect('equal')
        plt.gca().autoscale(tight=True)
        plt.colorbar()
        plt.show()

        plt.figure()
        plt.imshow(g.l)
        plt.title('l')
        plt.gca().set_aspect('equal')
        plt.gca().autoscale(tight=True)
        plt.colorbar()
        plt.show()

        plt.figure()
        plt.imshow(g.m)
        plt.title('m')
        plt.gca().set_aspect('equal')
        plt.gca().autoscale(tight=True)
        plt.colorbar()
        plt.show()

    if timeForLoop:
        dx = 50.0
        x = np.arange(0.0, 20000.1, dx)
        z = np.arange(0.0, 15000.1, dx)

        g = GrilleFDTD(x, z)

        Vp = 4000.0 + np.zeros((z.size, x.size))
        Vp[200:,:] = 5000.0
        sigma = 0.25              # coeff Poisson
        Vs = Vp * np.sqrt((0.5-sigma)/(1.0-sigma))
        rho = 2670.0 + np.zeros(Vp.shape)
        rho[200:,:] = 2700.0
        dt = 0.006

        g.defProp(Vp, Vs, rho, dt)

        src = Ricker(200, 150, 5, dt)

        t = 0.5

        tic = time.time()
        g.propage0(src, t, True)
        t0 = time.time() - tic
        print(t0)
        tic = time.time()
        g.propage(src, t, True)
        t = time.time() - tic
        print(t)
