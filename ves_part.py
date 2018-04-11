#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Programmes pour la modélisation et l'inversion de sondages géoélectriques
"""
import numpy as np
import matplotlib.pyplot as plt


class Sondage:
    """
    Classe pour modéliser et inverser les sondages géoélectriques en
    dispositif Schlumberger
    """
    def __init__(self, data):
        """
        Constructeur

        Parameters
        ----------
        data : données, array numpy N x 2
             1re colonne: AB/2 (moitié de la distance entre les électrodes
                                de courant)
             2e colonne : résistivité apparente mesurée

        """
        self.ab2 = data[:, 0]
        self.rhoa = data[:, 1]
        self.rhoam = np.ones(self.rhoa.shape)
        q = 13
        self.f = 10.0
        self.m = 4.438
        self.x = 0.0
        self.e = np.exp(0.5*np.log(10.0)/self.m)
        h = 2*q-2
        n = 1
        self.a = np.empty((n+h,))

    @property
    def data(self):
        return np.vstack((self.ab2, self.rhoa)).T

    @data.setter
    def data(self, data):
        self.ab2 = data[:, 0]
        self.rhoa = data[:, 1]
        self.rhoam = np.ones(self.rhoa.shape)

    def mod(self, m, rhoa=None):
        """
        Modélisation d'un sondage Schlumberger

        Parameters
        ----------
        m : modèle, [résistivités des couches, épaisseurs des couches]
              l'épaisseur de la dernière couche est infinie et est omise, e.g.
              m = np.array([rho_0, rho_1, rho_2, e_0, e_1])
        rhoa : array numpy préalloué (un array est créé si None)

        Returns
        -------
        rhoa : résistivités apparentes (array numpy)
        """
        if rhoa is None:
            rhoa = np.empty(self.rhoa.shape)
        lr = 1+int(m.size/2)
        r = m[:lr]
        t = m[lr:]
        for ns in np.arange(self.ab2.size):
            u = self.ab2[ns]*np.exp(-self.f*np.log(10)/self.m-self.x)
            for i in np.arange(self.a.size):
                w = r.size-1
                v = r[r.size-1]
                while w > 0:
                    w = w-1
                    aa = np.tanh(t[w]/u)
                    v = (v+r[w]*aa)/(1+v*aa/r[w])
                self.a[i] = v
                u = u*self.e

            rhoa[ns] = 105.0*self.a[0] - 262.0*self.a[2] + 416.0*self.a[4]
            rhoa[ns] += -746*self.a[6] + 1605*self.a[8]
            rhoa[ns] += -4390.0*self.a[10] + 13396.0*self.a[12]
            rhoa[ns] += -27841.0*self.a[14] + 16448.0*self.a[16]
            rhoa[ns] += 8183.0*self.a[18] + 2525.0*self.a[20]
            rhoa[ns] = (rhoa[ns] + 336.0*self.a[22] + 225.0*self.a[24])/10000.0
        return rhoa

    def _jacobian(self, m):
        par = 0.1
        m2 = m.copy()
        rhoa2 = np.empty((self.ab2.size, ))
        A = np.empty((self.ab2.size, m.size))
        for i2 in np.arange(m.size):
            m2[i2] = (m[i2]*par)+m[i2]
            rhoa2 = self.mod(m2)
            # les termes de la jacobienne sont normalisés
            A[:, i2] = ((rhoa2-self.rhoam)/(m[i2]*par))*m[i2]/self.rhoa
            m2[:] = m[:]

        return A

    def inv(self, m0, maxit=100, eps=1.e-8, showfig=False):
        """
        Inversion d'un sondage Schlumberger

        Parameters
        ----------
        m0 : modèle initial [résistivités des couches, épaisseurs des couches]
              l'épaisseur de la dernière couche est infinie et est omise, e.g.
              m0 = np.array([rho_0, rho_1, rho_2, e_0, e_1])
        maxit : nombre maximum d'itérations
        eps : critère de convergence, arrêt si |E^(p) - E^(p-1)| < eps*E^(p-1)
        showfig : affiche le modèle et l'ajustement durant l'inversion

        Returns
        -------
        m : modèle estimé (array numpy)
        """

        m = m0.copy()
        kr = eps
        iteration = 1
        maxiteration = maxit
        dfit = 1

        if showfig:
            fig, ax = plt.subplots(figsize=(8, 8), ncols=2)
            l1, = ax[0].loglog(self.ab2, self.rhoa, 'o')
            l2, = ax[0].loglog(self.ab2, self.rhoam)
            ax[0].set_xlabel('AB/2 (m)', fontsize=16)
            ax[0].set_ylabel('$\\rho_a (\Omega m)$', usetex=True, fontsize=16)
            lr = 1+int(m.size/2)
            z = np.cumsum(m[lr:])
            z = np.kron(z, np.ones((2,)))
            z = np.hstack((0.0, z, 2*z[-1]))
            rho = np.kron(m[:lr], np.ones((2,)))
            l3, = ax[1].semilogx(rho, z)
            ax[1].set_xlabel('$\\rho (\Omega m)$', usetex=True, fontsize=16)
            ax[1].set_ylabel('Profondeur (m)', fontsize=16)
            ax[1].invert_yaxis()
            ax[1].yaxis.set_label_position('right')
            ax[1].yaxis.tick_right()
            stitle = fig.suptitle('')
            plt.show(block=False)

        while iteration < maxiteration:
            self.mod(m, self.rhoam)

            e1 =
            e1 = e1.reshape(-1, 1)
            misfit1 = e1.T.dot(e1)
            if misfit1 < kr:
                return m
            J = self._jacobian(m)
            U, s, Vh = np.linalg.svd(J, full_matrices=False)
            V = Vh.T
            l = 0
            k = 1
            S = np.eye(s.size)
            while l < s.size:
                beta =
                if beta < 1.e-5:
                    beta = 0.001*(l+1)
                for i4 in np.arange(s.size):
                    S[i4, i4] =
                dmg =
                mg =
                self.mod(mg, self.rhoam)
                e2 =
                misfit2 = e2.T.dot(e2)
                if misfit2 > misfit1:
                    l += 1
                    k += 1
                    if k == s.size:
                        return m
                else:
                    l = s.size+1
                    m = mg
                    dfit = (misfit1-misfit2)/misfit1
                    iteration += 1
                    if showfig:
                        l2.set_ydata(self.rhoam)
                        z = np.cumsum(m[lr:])
                        z = np.kron(z, np.ones((2,)))
                        z = np.hstack((0.0, z, 2*z[-1]))
                        rho = np.kron(m[:lr], np.ones((2,)))
                        l3.set_xdata(rho)
                        l3.set_ydata(z)
                        ax[0].relim()
                        ax[0].autoscale()
                        ax[1].relim()
                        ax[1].autoscale()
                        stitle.set_text('Itération no {0:d}'.format(iteration))
                        fig.canvas.draw()
                        plt.pause(0.1)
                    if dfit < kr:
                        return m
        return m
