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

class CPML:
    def __init__(self, N=20, nd=2.0, Rc=0.001, nu0=2.0, nnu=2.0,
                 alpha0=20*np.pi, nalpha=1.0):
        """
        Input
            N      : nombre de couches PML
            nd     : ordre du profile du facteur d'amortissement
            Rc     : coefficient de réflexion théorique à la limite des PML
            nu0    : valeur max du paramètre nu
            nnu    : ordre du profile du paramètre nu
            nalpha : ordre du profile du paramètre alpha
            alpha0 : valeur max du paramètre alpha
        """
        self.N = N
        self.nd = nd
        self.Rc = Rc
        self.nu0 = nu0
        self.nnu = nnu
        self.alpha0 = alpha0
        self.nalpha = nalpha

    def prepare(self, nx, nz, dx, dt, V):
        """
        Calcul des coefficients utilisés avec les PML (taille: nz x nx)

          inux     : 1/nu_x évalué à x = i Delta x
          bx       : b_x    évalué à x = i Delta x
          cx       : c_x    évalué à x = i Delta x
          inux2    : 1/nu_x évalué à x = (i+1/2) Delta x
          bx2      : b_x    évalué à x = (i+1/2) Delta x
          cx2      : c_x    évalué à x = (i+1/2) Delta x
          inuz     : 1/nu_z évalué à z = j Delta z
          bz       : b_z    évalué à z = j Delta z
          cz       : c_z    évalué à z = j Delta z
          inuz2    : 1/nu_z évalué à z = (j+1/2) Delta z
          bz2      : b_z    évalué à z = (j+1/2) Delta z
          cz2      : c_z    évalué à z = (j+1/2) Delta z
          psi_txxx : variable mémoire pour d tau_xx / dx
          psi_txzz : variable mémoire pour d tau_xz / dz
          psi_txzx : variable mémoire pour d tau_xz / dx
          psi_tzzz : variable mémoire pour d tau_zz / dz
          psi_vxx  : variable mémoire pour d v_x / dx
          psi_vzz  : variable mémoire pour d v_z / dz
          psi_vxz  : variable mémoire pour d v_x / dz
          psi_vzx  : variable mémoire pour d v_z / dx

        """
        # à (i,j)

        xp =
        d0 = (self.nd+1) * np.log(1/self.Rc)*V / (2*self.N*dx)
        dx_pml =

        nu = 1. + (self.nu0-1.) * (xp / (self.N*dx))**self.nnu
        alpha_pml = self.alpha0*(1. - (xp / (self.N*dx))**self.nalpha)

        self.inux = np.ones((nz, nx))
        self.inux[] =
        self.inux[] =

        d = np.zeros((nz, nx))
        d[] =
        d[] =

        alpha = np.zeros((nz, nx))
        alpha[] =
        alpha[] =

        self.bx = np.exp(-(d * self.inux + alpha)*dt)
        with np.errstate(divide='ignore', invalid='ignore'):
            self.cx = d * self.inux * (self.bx-1.) / (d+alpha / self.inux)
        self.cx[np.isnan(self.cx)] = 0.0

        self.inuz = np.ones((nz, nx))
        self.inuz[] =
        self.inuz[] =

        d = np.zeros((nz, nx))
        d[] =
        d[] =

        alpha = np.zeros((nz, nx))
        alpha[] =
        alpha[] =

        self.bz = np.exp(-(d * self.inuz + alpha)*dt)
        with np.errstate(divide='ignore', invalid='ignore'):
            self.cz = d * self.inuz * (self.bz-1.) / (d+alpha / self.inuz)
        self.cz[np.isnan(self.cz)] = 0.0

        # à (i+1/2,j+1/2)

        xp1 =
        xp2 =
        dx1_pml = d0 * (xp1 / (self.N*dx))**self.nd
        dx2_pml = d0 * (xp2 / (self.N*dx))**self.nd

        nu1 = 1. + (self.nu0-1) * (xp1 / (self.N*dx))**self.nnu
        nu2 = 1. + (self.nu0-1) * (xp2 / (self.N*dx))**self.nnu
        alpha1_pml = self.alpha0 * (1. - (xp1 / (self.N*dx))**self.nalpha)
        alpha2_pml = self.alpha0 * (1. - (xp2 / (self.N*dx))**self.nalpha)

        self.inux2 = np.ones((nz, nx))
        self.inux2[] =
        self.inux2[] =

        d = np.zeros((nz, nx))
        d[] =
        d[] =

        alpha = np.zeros((nz, nx))
        alpha[] =
        alpha[] =

        self.bx2 = np.exp(-(d * self.inux2 + alpha)*dt)
        with np.errstate(divide='ignore', invalid='ignore'):
            self.cx2 = d * self.inux2 * (self.bx2-1.) / (d+alpha / self.inux2)
        self.cx2[np.isnan(self.cx2)] = 0.0

        self.inuz2 = np.ones((nz, nx))
        self.inuz2[] =
        self.inuz2[] =

        d = np.zeros((nz, nx))
        d[] =
        d[] =

        alpha = np.zeros((nz, nx))
        alpha[] =
        alpha[] =

        self.bz2 = np.exp(-(d * self.inuz2 + alpha)*dt)
        with np.errstate(divide='ignore', invalid='ignore'):
            self.cz2 = d * self.inuz2 * (self.bz2-1.) / (d+alpha / self.inuz2)
        self.cz2[np.isnan(self.cz2)] = 0.0

        self.psi_txxx = np.zeros((nz, nx))
        self.psi_txzz = np.zeros((nz, nx))
        self.psi_txzx = np.zeros((nz, nx))
        self.psi_tzzz = np.zeros((nz, nx))
        self.psi_vxx = np.zeros((nz, nx))
        self.psi_vzz = np.zeros((nz, nx))
        self.psi_vxz = np.zeros((nz, nx))
        self.psi_vzx = np.zeros((nz, nx))

    def reset_vm(self):
        self.psi_txxx[:,:] = 0.0
        self.psi_txzz[:,:] = 0.0
        self.psi_txzx[:,:] = 0.0
        self.psi_tzzz[:,:] = 0.0
        self.psi_vxx[:,:] = 0.0
        self.psi_vzz[:,:] = 0.0
        self.psi_vxz[:,:] = 0.0
        self.psi_vzx[:,:] = 0.0


class A1():
    def __init__(self, Vp, Vs, dx, dt):
        self.alpha_h = Vp[0, :] * dt/dx
        self.beta_h  = Vs[0, :] * dt/dx
        self.alpha_b = Vp[-1,:] * dt/dx
        self.beta_b  = Vs[-1,:] * dt/dx
        self.alpha_g = Vp[:, 0] * dt/dx
        self.beta_g  = Vs[:, 0] * dt/dx
        self.alpha_d = Vp[:,-1] * dt/dx
        self.beta_d  = Vs[:,-1] * dt/dx


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
                ax1.clear()
                ax1.imshow(v_x)
                ax2.clear()
                ax2.imshow(v_z)
                fig.canvas.draw()
                stitle.set_text('t = {0:6.3f} (it no {1:d}/{2:d})'.\
                                format(m*self.dt, m+1, nstep))
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
                ax1.clear()
                ax1.imshow(v_x)
                ax2.clear()
                ax2.imshow(v_z)
                fig.canvas.draw()
                stitle.set_text('t = {0:6.3f} (it no {1:d}/{2:d})'.\
                                format(m*self.dt, m+1, nstep))
                plt.pause(0.01)

    def propageO24(self, src, t, showPlot=False):
        nstep = 1 + int(t/self.dt)

        c1 = 9./8.
        c2 = 1./24.

        tau_xx = np.zeros((self.nz, self.nx))
        tau_zz = np.zeros((self.nz, self.nx))
        tau_xz = np.zeros((self.nz, self.nx))
        v_x = np.zeros((self.nz, self.nx))
        v_z = np.zeros((self.nz, self.nx))

        if showPlot:
            fig, (ax1, ax2) = plt.subplots(ncols=2)
            stitle = fig.suptitle('')
            plt.show(block=False)

        for m in range(nstep):

            # applique la source au noeud (i+1/2,j)
            tau_zz[src.j, src.i] += src(m)
            tau_xx[src.j, src.i] += src(m)

            v_x[] += self.b1[] * (
                    c1 * (tau_xx[] - tau_xx[])-
                    c2 * (tau_xx[] - tau_xx[])+
                    c1 * (tau_xz[] - tau_xz[])-
                    c2 * (tau_xz[] - tau_xz[]))

            v_z[] += self.b2[] * (
                    c1 * (tau_xz[] - tau_xz[])-
                    c2 * (tau_xz[] - tau_xz[])+
                    c1 * (tau_zz[] - tau_zz[])-
                    c2 * (tau_zz[] - tau_zz[]))

            tau_xx[] += self.lm[] * (
                    c1 * (v_x[] - v_x[]) -
                    c2 * (v_x[] - v_x[])) + \
                    self.l[] * (
                    c1 * (v_z[] - v_z[]) -
                    c2 * (v_z[] - v_z[]))

            tau_zz[] += self.lm[] * (
                    c1 * (v_z[] - v_z[]) -
                    c2 * (v_z[] - v_z[])) + \
                    self.l[] * (
                    c1 * (v_x[] - v_x[]) -
                    c2 * (v_x[] - v_x[]))

            tau_xz[] += self.m[] * (
                    c1 * (v_x[] - v_x[]) -
                    c2 * (v_x[] - v_x[]) +
                    c1 * (v_z[] - v_z[]) -
                    c2 * (v_z[] - v_z[]))

            if showPlot and np.remainder(m, 20) == 0:
                ax1.clear()
                ax1.imshow(v_x)
                ax2.clear()
                ax2.imshow(v_z)
                fig.canvas.draw()
                stitle.set_text('t = {0:6.3f} (it no {1:d}/{2:d})'.\
                                format(m*self.dt, m+1, nstep))
                plt.pause(0.01)

    def propageO24_pml(self, src, t, pml=None, a1=None, dirichlet=False,
                       trc=[], calcE=False, showPlot=False):
        nstep = 1 + int(t/self.dt)

        if calcE is True:
            E = np.zeros((nstep,))
            rho1 = self.dt / (self.b1*self.dx)
            rho2 = self.dt / (self.b2*self.dx)
            lambda2mu = self.lm*self.dx/self.dt
            lambd = self.l*self.dx/self.dt
            mu = self.m*self.dx/self.dt

            ind = np.ones(mu.shape, dtype=np.bool)
            if pml is not None:
                ind[:pml.N,:] = False
                ind[:,:pml.N] = False
                ind[-pml.N:,:] = False
                ind[:,-pml.N:] = False

        if len(trc) > 0:
            ntrc = trc.shape[0]
            trc_vx = np.zeros((nstep, ntrc))
            trc_vz = np.zeros((nstep, ntrc))
            # on arrondi à l'indice le plus proche
            itrc = np.array(np.round((trc[:,0]-self.x[0])/self.dx),
                            dtype=np.int64)
            jtrc = np.array(np.round((trc[:,1]-self.z[0])/self.dx),
                            dtype=np.int64)

        c1 = 9./8.
        c2 = 1./24.

        tau_xx = np.zeros((self.nz, self.nx))
        tau_zz = np.zeros((self.nz, self.nx))
        tau_xz = np.zeros((self.nz, self.nx))
        v_x = np.zeros((self.nz, self.nx))
        v_z = np.zeros((self.nz, self.nx))

        if pml is not None:
            b1 = self.b1 * self.dx
            b2 = self.b2 * self.dx
            lm = self.lm * self.dx
            ll = self.l  * self.dx
            mm = self.m  * self.dx

            c1 /= self.dx
            c2 /= self.dx
        else:
            b1 = self.b1
            b2 = self.b2
            lm = self.lm
            ll = self.l
            mm = self.m

        if showPlot:
            fig, (ax1, ax2) = plt.subplots(ncols=2)
            stitle = fig.suptitle('')
            plt.show(block=False)

        for m in range(nstep):

            # applique la source au noeud (i+1/2,j)
            tau_zz[src.j, src.i] += src(m)
            tau_xx[src.j, src.i] += src(m)

            if pml is not None and a1 is not None:
                # A1

                # en haut
                v_x[1,2:-2]=(1-a1.beta_h[2:-2]) * v_x[1,2:-2]+a1.beta_h[2:-2] * v_x[2,2:-2]
                v_x[0,2:-2]=(1-a1.beta_h[2:-2]) * v_x[0,2:-2]+a1.beta_h[2:-2] * v_x[1,2:-2]
                # en bas
                v_x[-2,2:-2]=(1-a1.beta_b[2:-2]) * v_x[-2,2:-2]+a1.beta_b[2:-2] * v_x[-3,2:-2]
                v_x[-1,2:-2]=(1-a1.beta_b[2:-2]) * v_x[-1,2:-2]+a1.beta_b[2:-2] * v_x[-2,2:-2]
                # à droite
                v_x[2:-2,-2]=(1-a1.alpha_d[2:-2]) * v_x[2:-2,-2]+a1.alpha_d[2:-2] * v_x[2:-2,-3]
                v_x[2:-2,-1]=(1-a1.alpha_d[2:-2]) * v_x[2:-2,-1]+a1.alpha_d[2:-2] * v_x[2:-2,-2]
                # à gauche
                v_x[2:-2,1]=(1-a1.alpha_g[2:-2]) * v_x[2:-2,1]+a1.alpha_g[2:-2] * v_x[2:-2,2]
                v_x[2:-2,0]=(1-a1.alpha_g[2:-2]) * v_x[2:-2,0]+a1.alpha_g[2:-2] * v_x[2:-2,1]

                # PML -> update nodes that were not updated with A1
                ddx=c1*(tau_xx[2:-2,2:-2]-tau_xx[2:-2,1:-3])-c2*(tau_xx[2:-2,3:-1]-tau_xx[2:-2,:-4])
                ddz=c1*(tau_xz[2:-2,2:-2]-tau_xz[1:-3,2:-2])-c2*(tau_xz[3:-1,2:-2]-tau_xz[:-4,2:-2])

                pml.psi_txxx[2:-2,2:-2] = pml.bx[2:-2,2:-2] * pml.psi_txxx[2:-2,2:-2] + pml.cx[2:-2,2:-2] * ddx
                pml.psi_txzz[2:-2,2:-2] = pml.bz[2:-2,2:-2] * pml.psi_txzz[2:-2,2:-2] + pml.cz[2:-2,2:-2] * ddz

                ddx=pml.inux[2:-2,2:-2] * ddx + pml.psi_txxx[2:-2,2:-2]
                ddz=pml.inuz[2:-2,2:-2] * ddz + pml.psi_txzz[2:-2,2:-2]

                v_x[2:-2,2:-2] += b1[2:-2,2:-2] * (ddx + ddz)

            elif pml is not None:
                ddx=c1*(tau_xx[]-tau_xx[])-c2*(tau_xx[]-tau_xx[])
                ddz=c1*(tau_xz[]-tau_xz[])-c2*(tau_xz[]-tau_xz[])

                pml.psi_txxx[] = pml.bx[]*pml.psi_txxx[] + pml.cx[]*ddx
                pml.psi_txzz[] = pml.bz[]*pml.psi_txzz[] + pml.cz[]*ddz

                ddx=pml.inux[]*ddx + pml.psi_txxx[]
                ddz=pml.inuz[]*ddz + pml.psi_txzz[]

                v_x[] += b1[]*(ddx + ddz)

            elif a1 is not None:
                # en haut
                v_x[1,2:-2]=(1-a1.beta_h[2:-2]) * v_x[1,2:-2]+a1.beta_h[2:-2] * v_x[2,2:-2]
                v_x[0,2:-2]=(1-a1.beta_h[2:-2]) * v_x[0,2:-2]+a1.beta_h[2:-2] * v_x[1,2:-2]
                # en bas
                v_x[-2,2:-2]=(1-a1.beta_b[2:-2]) * v_x[-2,2:-2]+a1.beta_b[2:-2] * v_x[-3,2:-2]
                v_x[-1,2:-2]=(1-a1.beta_b[2:-2]) * v_x[-1,2:-2]+a1.beta_b[2:-2] * v_x[-2,2:-2]
                # à droite
                v_x[2:-2,-2]=(1-a1.alpha_d[2:-2]) * v_x[2:-2,-2]+a1.alpha_d[2:-2] * v_x[2:-2,-3]
                v_x[2:-2,-1]=(1-a1.alpha_d[2:-2]) * v_x[2:-2,-1]+a1.alpha_d[2:-2] * v_x[2:-2,-2]
                # à gauche
                v_x[2:-2,1]=(1-a1.alpha_g[2:-2]) * v_x[2:-2,1]+a1.alpha_g[2:-2] * v_x[2:-2,2]
                v_x[2:-2,0]=(1-a1.alpha_g[2:-2]) * v_x[2:-2,0]+a1.alpha_g[2:-2] * v_x[2:-2,1]

                v_x[2:-2,2:-2] += b1[2:-2,2:-2] * (
                    c1*(tau_xx[2:-2,2:-2]-tau_xx[2:-2,1:-3])-c2*(tau_xx[2:-2,3:-1]-tau_xx[2:-2,:-4])+
                    c1*(tau_xz[2:-2,2:-2]-tau_xz[1:-3,2:-2])-c2*(tau_xz[3:-1,2:-2]-tau_xz[:-4,2:-2]))

            else:
                v_x[] += b1[] * (
                        c1 * (tau_xx[]-tau_xx[])-
                        c2 * (tau_xx[] -tau_xx[])+
                        c1 * (tau_xz[]-tau_xz[])-
                        c2 * (tau_xz[] -tau_xz[]))

            # Vz
            if pml is not None and a1 is not None:
                # A1
                # en haut
                v_z[1,2:-2]=(1-a1.alpha_h[2:-2]) * v_z[1,2:-2]+a1.alpha_h[2:-2] * v_z[2,2:-2]
                v_z[0,2:-2]=(1-a1.alpha_h[2:-2]) * v_z[0,2:-2]+a1.alpha_h[2:-2] * v_z[1,2:-2]
                # en bas
                v_z[-2,2:-2]=(1-a1.alpha_b[2:-2]) * v_z[-2,2:-2]+a1.alpha_b[2:-2] * v_z[-3,2:-2]
                v_z[-1,2:-2]=(1-a1.alpha_b[2:-2]) * v_z[-1,2:-2]+a1.alpha_b[2:-2] * v_z[-2,2:-2]
                # à droite
                v_z[2:-2,-2]=(1-a1.beta_d[2:-2]) * v_z[2:-2,-2]+a1.beta_d[2:-2] * v_z[2:-2,-3]
                v_z[2:-2,-1]=(1-a1.beta_d[2:-2]) * v_z[2:-2,-1]+a1.beta_d[2:-2] * v_z[2:-2,-2]
                # à gauche
                v_z[2:-2,1]=(1-a1.beta_d[2:-2]) * v_z[2:-2,1]+a1.beta_d[2:-2] * v_z[2:-2,2]
                v_z[2:-2,0]=(1-a1.beta_d[2:-2]) * v_z[2:-2,0]+a1.beta_d[2:-2] * v_z[2:-2,1]

                # PML
                ddx=c1*(tau_xz[2:-2,3:-1]-tau_xz[2:-2,2:-2])-c2*(tau_xz[2:-2,4:]-tau_xz[2:-2,1:-3])
                ddz=c1*(tau_zz[3:-1,2:-2]-tau_zz[2:-2,2:-2])-c2*(tau_zz[4:,2:-2]-tau_zz[1:-3,2:-2])

                pml.psi_txzx[2:-2,2:-2] = pml.bx2[2:-2,2:-2] * pml.psi_txzx[2:-2,2:-2] + pml.cx2[2:-2,2:-2] * ddx
                pml.psi_tzzz[2:-2,2:-2] = pml.bz2[2:-2,2:-2] * pml.psi_tzzz[2:-2,2:-2] + pml.cz2[2:-2,2:-2] * ddz

                ddx=pml.inux2[2:-2,2:-2] * ddx + pml.psi_txzx[2:-2,2:-2]
                ddz=pml.inuz2[2:-2,2:-2] * ddz + pml.psi_tzzz[2:-2,2:-2]

                v_z[2:-2,2:-2] += b2[2:-2,2:-2] * (ddx + ddz)

            elif pml is not None:
                ddx=c1*(tau_xz[]-tau_xz[])-c2*(tau_xz[]-tau_xz[1:-2,:-3])
                ddz=c1*(tau_zz[]-tau_zz[])-c2*(tau_zz[3:,1:-2]-tau_zz[])

                pml.psi_txzx[] = pml.bx2[]*pml.psi_txzx[] + pml.cx2[]*ddx

                pml.psi_tzzz[] = pml.bz2[]*pml.psi_tzzz[] + pml.cz2[]*ddz

                ddx=pml.inux2[]*ddx + pml.psi_txzx[]
                ddz=pml.inuz2[]*ddz + pml.psi_tzzz[]

                v_z[] += b2[]*(ddx + ddz)

            elif a1 is not None:
                # en haut
                v_z[1,2:-2]=(1-a1.alpha_h[2:-2]) * v_z[1,2:-2]+a1.alpha_h[2:-2] * v_z[2,2:-2]
                v_z[0,2:-2]=(1-a1.alpha_h[2:-2]) * v_z[0,2:-2]+a1.alpha_h[2:-2] * v_z[1,2:-2]
                # en bas
                v_z[-2,2:-2]=(1-a1.alpha_b[2:-2]) * v_z[-2,2:-2]+a1.alpha_b[2:-2] * v_z[-3,2:-2]
                v_z[-1,2:-2]=(1-a1.alpha_b[2:-2]) * v_z[-1,2:-2]+a1.alpha_b[2:-2] * v_z[-2,2:-2]
                # à droite
                v_z[2:-2,-2]=(1-a1.beta_d[2:-2]) * v_z[2:-2,-2]+a1.beta_d[2:-2] * v_z[2:-2,-3]
                v_z[2:-2,-1]=(1-a1.beta_d[2:-2]) * v_z[2:-2,-1]+a1.beta_d[2:-2] * v_z[2:-2,-2]
                # à gauche
                v_z[2:-2,1]=(1-a1.beta_d[2:-2]) * v_z[2:-2,1]+a1.beta_d[2:-2] * v_z[2:-2,2]
                v_z[2:-2,0]=(1-a1.beta_d[2:-2]) * v_z[2:-2,0]+a1.beta_d[2:-2] * v_z[2:-2,1]

                v_z[2:-2,2:-2] += b2[2:-2,2:-2] * (
                    c1*(tau_xz[2:-2,3:-1]-tau_xz[2:-2,2:-2])-c2*(tau_xz[2:-2,4:]-tau_xz[2:-2,1:-3])+
                    c1*(tau_zz[3:-1,2:-2]-tau_zz[2:-2,2:-2])-c2*(tau_zz[4:,2:-2]-tau_zz[1:-3,2:-2]))

            else:
                v_z[] += b2[] * (
                        c1 * (tau_xz[] - tau_xz[])-
                        c2 * (tau_xz[] - tau_xz[1:-2,:-3])+
                        c1 * (tau_zz[] - tau_zz[])-
                        c2 * (tau_zz[3:,1:-2] - tau_zz[]))

            if len(trc) > 0:
                for nt in np.arange(ntrc):
                    trc_vx[m, nt] = v_x[jtrc[nt], itrc[nt]]
                    trc_vz[m, nt] = v_z[jtrc[nt], itrc[nt]]

            # source alternative
            #v_z[src.j, src.i] += dx * b2[src.j, src.i] * src(m)
            # v_x[src.j, src.i] += dx * b1[src.j, src.i] * src(m)

            if dirichlet is True:
                # Conditions de Dirichlet: on prend 2 rangées à cause de l'ordre 4 de
                # l'opérateur spatial
                v_x[:,:2]  = 0.0
                v_x[:,-2:] = 0.0
                v_x[:2,:]  = 0.0
                v_x[-2:,:] = 0.0
                v_z[:,:2]  = 0.0
                v_z[:,-2:] = 0.0
                v_z[:2,:]  = 0.0
                v_z[-2:,:] = 0.0

            # Contraintes
            if pml is not None:
                ddx=c1*(v_x[]-v_x[])-c2*(v_x[]-v_x[])
                ddz=c1*(v_z[]-v_z[])-c2*(v_z[3:,1:-2]-v_z[])

                pml.psi_vxx[] = pml.bx2[]*pml.psi_vxx[] + pml.cx2[]*ddx

                pml.psi_vzz[] = pml.bz[]*pml.psi_vzz[] + pml.cz[]*ddz

                ddx=pml.inux2[]*ddx + pml.psi_vxx[]
                ddz=pml.inuz[]*ddz + pml.psi_vzz[]

                tau_xx[] += lm[]*ddx + ll[]*ddz

                tau_zz[] += lm[]*ddz + ll[]*ddx

                ddx=c1*(v_z[]-v_z[])-c2*(v_z[]-v_z[1:-2,:-3])
                ddz=c1*(v_x[]-v_x[])-c2*(v_x[]-v_x[])

                pml.psi_vzx[] = pml.bx[]*pml.psi_vzx[] + pml.cx[]*ddx

                pml.psi_vxz[] = pml.bz2[]*pml.psi_vxz[] + pml.cz2[]*ddz

                ddx=pml.inux[]*ddx + pml.psi_vzx[]
                ddz=pml.inuz2[]*ddz + pml.psi_vxz[]

                tau_xz[] += mm[]*(ddx + ddz)

            else:
                tau_xx[] += lm[] * (
                        c1 * (v_x[] - v_x[]) -
                        c2 * (v_x[] - v_x[])) + \
                        ll[] * (
                        c1 * (v_z[] - v_z[]) -
                        c2 * (v_z[3:,1:-2] - v_z[]))

                tau_zz[] += lm[] * (
                        c1 * (v_z[] - v_z[]) -
                        c2 * (v_z[3:,1:-2] - v_z[])) + \
                        ll[] * (
                        c1 * (v_x[] - v_x[]) -
                        c2 * (v_x[] - v_x[]))

                tau_xz[] += mm[] * (
                        c1 * (v_x[] - v_x[]) -
                        c2 * (v_x[] - v_x[]) +
                        c1 * (v_z[] - v_z[]) -
                        c2 * (v_z[] - v_z[1:-2,:-3]))

            if calcE:
                tmp = 0.5 * (rho1 * v_x**2 + rho2 * v_z**2)
                epsilon_xx = (lambda2mu * tau_xx - lambd * tau_zz) / (4*mu*(lambd + mu))
                epsilon_yy = (lambda2mu * tau_zz - lambd * tau_xx) / (4*mu*(lambd + mu))
                epsilon_xy = tau_xz / (2 * mu)
                tmp += 0.5 * (epsilon_xx*tau_xx + epsilon_yy*tau_zz + 2*epsilon_xy*tau_xz)
                E[m] = np.sum(tmp[ind])

            if showPlot and np.remainder(m, 20) == 0:
                ax1.clear()
                ax1.imshow(v_x)
                ax2.clear()
                ax2.imshow(v_z)
                fig.canvas.draw()
                stitle.set_text('t = {0:6.3f} (it no {1:d}/{2:d})'.format(m*self.dt, m+1, nstep))
                plt.pause(0.01)

        if calcE and len(trc) > 0:
            return trc_vx, trc_vz, E
        elif calcE:
            return E
        elif len(trc) > 0:
            return trc_vx, trc_vz


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
