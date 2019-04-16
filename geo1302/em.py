#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Codes de modélisation EM

Ces codes sont une adaptation des codes matlab disponibles à

https://github.com/brianborchers/PEIP

Created on Mon Apr  8 19:11:11 2019

@author: giroux
"""
import os
import numpy as np
import matplotlib.pyplot as plt

dir_path = os.path.dirname(os.path.abspath(__file__))

MU0 = np.pi * 4.e-7
WT0 = np.load(dir_path+'/wt0.npy').flatten()
WT1 = np.load(dir_path+'/wt1.npy').flatten()
YBASE = np.load(dir_path+'/ybase.npy').flatten()


class Ground:
    """
    Classe pour décrire un sol stratifié

    Paramètres
    ----------
    nlayer : int
        nombre de couches
    d : float
        épaisseur des couches
    """
    def __init__(self, nlayer=11, d=0.2):
        self.nlayer = nlayer
        self.d = d * np.ones((nlayer-1,))
        self.mu = MU0 + np.zeros((nlayer,))


    def plot(self, sigma, ax=None):
        """
        Trace la conductivité en fct de la profondeur

        Paramètres
        ----------
        sigma : array numpy
            conductivité des couches en S/m
        """
        z = np.cumsum(self.d)
        z = np.kron(z, np.ones((2,)))
        z = np.hstack((0.0, z, z[-1]+z[0]))
        s = 1000*np.kron(sigma, np.ones((2,)))

        if ax is None:
            plt.figure(figsize=(6,8))
            ax = plt.gca()
        ax.plot(s, z)
        ax.invert_yaxis()
        ax.set_ylabel('Profondeur (m)', fontsize=14)
        ax.set_xlabel('Conductivité (mS/m)', fontsize=14)


class EM38:
    """
    Classe pour modéliser la réponse d'un conductivimètre EM-38

    Paramètres
    ----------
    h : array numpy
        hauteurs des mesures (m)
    ground : object ground
        description du sol
    f : float
        Fréquence d'opération (Hz)
    r : float
        distance entre les bobines
    """
    def __init__(self, h, ground, f=14600.0, r=1.0):
        self.h = h   # hauteur des mesures, m
        self.f = f   # fréquence, Hz
        self.omega = 2.0 * np.pi * self.f
        self.r = r   # distance entre les bobines, m
        self.sigma0 = 0.0
        self.g = ground


    def G(self, sigma):
        """
        Modèle direct

        Paramètres
        ----------
        sigma : array numpy
            conductivité des couches (S/m)
        """
        delta = np.sqrt(2/(sigma[0] * MU0 * self.omega))

        pred = np.empty((self.h.size*2,))
        for i in range(self.h.size):
            pred[i], pred[i+self.h.size] = self._predict(self.h[i], delta, sigma)

        return pred


    def jac(self, sigma):
        """
        Calcul de la jacobienne

        Paramètres
        ----------
        sigma : array numpy
            conductivité des couches (S/m)
        """
        J = np.zeros((self.h.size*2, self.g.nlayer))
        I = np.eye(self.g.nlayer)
        h = 1.0e-8
        tmp = self.G(sigma)
        for i in range(self.g.nlayer):
            J[:, i] = (self.G(sigma+h*I[:, i]) - tmp) / h
        return J


    def _predict(self, h, delta, sigma):
        g = YBASE / (self.r/delta)
        r0 = self._r0(g/delta, sigma)
        f0 = -r0 * g * g * np.exp(-2. * g * h / delta)
        f1 = -r0 * g * np.exp(-2. * g * h / delta)
        T0 = WT0.dot(f0) / (self.r/delta)
        T2 = WT1.dot(f1) / (self.r/delta)
        predv = np.imag(1 + T0 * (self.r/delta)**3) * 1000*4/(MU0 * self.omega * self.r * self.r)
        predh = np.imag(1 + T2 * (self.r/delta)**2) * 1000*4/(MU0 * self.omega * self.r * self.r)

        return predv, predh


    def _r0(self, lmbda, sigma):
        Y = np.empty((self.g.nlayer-1), dtype=complex)
        r0 = np.empty((lmbda.size,), dtype=complex)

        imo = 1j * self.g.mu * self.omega
        ismo = sigma * imo

        for j in np.arange(lmbda.size):
            u = np.sqrt(lmbda[j]**2 + ismo)
            N = u / imo

            tanhud = np.tanh(u[self.g.nlayer-2]*self.g.d[self.g.nlayer-2])
            Y[self.g.nlayer-2] = N[self.g.nlayer-2] * (N[self.g.nlayer-1] + N[self.g.nlayer-2] * tanhud) / \
                                                  (N[self.g.nlayer-2] + N[self.g.nlayer-1] * tanhud)
            for k in np.arange(self.g.nlayer-3, -1, -1):
                tanhud = np.tanh(u[k]*self.g.d[k])
                Y[k] = N[k] * (Y[k+1] + N[k] * tanhud) / (N[k] + Y[k+1] * tanhud)

            N0 = np.sqrt(lmbda[j]**2 + 1j * self.sigma0 * MU0 * self.omega) / (1j * MU0 * self.omega)
            r0[j] = (N0-Y[0]) / (N0+Y[0])
        return r0


if __name__ == '__main__':

    heights = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.50])
    mtrue = np.array([100, 95, 90, 95, 100, 130, 160, 200, 250, 310, 360])/1000.0

    ground = Ground()
    em38 = EM38(heights, ground)

    import time

    tic = time.time()
    pred = em38.G(mtrue)
    toc = time.time() - tic
    print(toc)

    tic = time.time()
    J = em38.jac(mtrue)
    toc = time.time() - tic
    print(toc)

    ground.plot(mtrue)
    plt.show()
