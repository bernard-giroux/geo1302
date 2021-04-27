# -*- coding: utf-8 -*-
"""
Codes de modélisation EM

Ces codes sont une adaptation des codes matlab disponibles à

https://github.com/brianborchers/PEIP

Ce module est une version cython du module em.py.  Le coeur de la fonction de
modélisation comporte deux boucles "for" imbriquées, avec un appel de fonction
de tangente hyperbolique, qui sont lentes à exécuter.  L'objectif de la version
cython est donc d'améliorer la performance.

Pour compiler ce code, la commande à exécuter est

python setup.py build_ext --inplace

Documentation cython pertinente pour comprendre ce code:

https://cython.readthedocs.io/en/latest/src/userguide/numpy_tutorial.html

@author: giroux
"""

import os
import numpy as np
import matplotlib.pyplot as plt

cimport numpy as np
cimport cython

cdef extern from "<cmath>":
    double complex tanh(double complex)
    double exp(double)
    double sqrt(double)


DTYPE = np.double
ctypedef np.double_t DTYPE_t

dir_path = os.path.dirname(os.path.abspath(__file__))
cdef double MU0 = np.pi * 4.e-7
WT0 = np.load(dir_path+'/wt0.npy').flatten()
WT1 = np.load(dir_path+'/wt1.npy').flatten()
YBASE = np.load(dir_path+'/ybase.npy').flatten()


cdef class Ground:
    '''
    Classe pour décrire un sol stratifié.

    Paramètres
    ----------
    nlayer : int
        nombre de couches
    d : float
        épaisseur des couches
    '''
    cdef int _nlayer
    cdef object _d
    cdef object _mu

    def __cinit__(self, nlayer=11, d=0.2):
        self._nlayer = nlayer
        self._d = d * np.ones((nlayer-1,))
        self._mu = MU0 + np.zeros((nlayer,))

    @property
    def nlayer(self):
        return self._nlayer

    @property
    def d(self):
        return self._d

    @property
    def mu(self):
        return self._mu


    def plot(self, sigma, ax=None, *args, **kwargs):
        """
        Trace la conductivité en fct de la profondeur.

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
            plt.figure(figsize=(6, 8))
            ax = plt.gca()
        lines = ax.plot(s, z, *args, **kwargs)
        if not ax.yaxis_inverted():
            ax.invert_yaxis()
        ax.set_ylabel('Profondeur (m)', fontsize=14)
        ax.set_xlabel('Conductivité (mS/m)', fontsize=14)
        return lines


cdef class EM38:
    '''
    Classe pour modéliser la réponse d'un conductivimètre EM-38.

    Paramètres
    ----------
    h : array numpy
        hauteurs des mesures (m)
    ground : object Ground
        description du sol
    f : float
        Fréquence d'opération (Hz)
    r : float
        distance entre les bobines
    '''
    cdef double[:] h
    cdef double f
    cdef double omega
    cdef double r
    cdef double sigma0
    cdef object g

    def __cinit__(self, double[:] h, ground, f=14600.0, r=1.0):
        self.h = h   # hauteur des mesures, m
        self.f = f   # fréquence, Hz
        self.omega = 2.0 * np.pi * self.f
        self.r = r   # distance entre les bobines, m
        self.sigma0 = 0.0
        self.g = ground


    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def G(self, double[:] sigma):
        """
        Modèle direct.

        Paramètres
        ----------
        sigma : array numpy
            conductivité des couches (S/m)

        Retourne
        --------
        pred : array numpy
            réponse prédite.  La 1re moitié des valeurs est la réponse pour des
            bobines verticales aux différentes hauteurs de mesure, et la 2e
            moitié est la réponse pour des bobines horizontales.
        """
        cdef double delta = sqrt(2/(sigma[0] * MU0 * self.omega))
        _pred = np.empty((self.h.size*2,))
        cdef double[::1] pred = _pred
        cdef int n = self.h.size
        cdef int i
        for i in range(n):
            pred[i], pred[i+n] = self._predict(self.h[i], delta, sigma)

        return _pred


    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def jac(self, double[:] sigma):
        """
        Calcul de la jacobienne.

        Paramètres
        ----------
        sigma : array numpy
            conductivité des couches (S/m)

        Retourne
        --------
        J : array numpy
            matrice jacobienne
        """
        _J = np.zeros((self.h.size*2, self.g.nlayer))
        cdef double[:, ::1] J = _J
        # _I = h*np.eye(self.g.nlayer)
        # cdef double[:, :] I = _I
        cdef double h = 1.0e-8
        cdef double[::1] tmp = self.G(sigma)
        cdef int nlayer = self.g.nlayer
        cdef int i, ii
        _tmp2 = np.empty((sigma.size,))
        cdef double[::1] tmp2 = _tmp2
        for i in range(nlayer):
            for ii in range(sigma.size):
                tmp2[ii] = sigma[ii]
            tmp2[i] = sigma[i] + h
            tmp2 = self.G(tmp2)
            for ii in range(sigma.size):
                J[ii, i] = (tmp2[ii] - tmp[ii]) / h
        return _J


    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def _predict(self, double h, double delta, double[:] sigma):
        cdef int i
        _g = YBASE / (self.r/delta)
        cdef double[::1] g = _g
        _tmp = np.empty((g.size,))
        cdef double[::1] tmp = _tmp
        for i in range(g.size):
            tmp[i] = g[i]/delta
        cdef double complex[::1] r0 = self._r0(tmp, sigma)
        _f0 = np.empty((r0.size,), dtype=complex)
        _f1 = np.empty((r0.size,), dtype=complex)
        cdef double complex[::1] f0 = _f0
        cdef double complex[::1] f1 = _f1
        for i in range(g.size):
            f0[i] = -r0[i] * g[i] * g[i] * exp(-2. * g[i] * h / delta)
            f1[i] = -r0[i] * g[i] * exp(-2. * g[i] * h / delta)
        cdef double complex T0 = 0.0
        cdef double complex T2 = 0.0
        cdef double[::1] wt0 = WT0
        cdef double[::1] wt1 = WT1
        for i in range(f0.size):
            T0 += wt0[i] * f0[i]
            T2 += wt1[i] * f1[i]
        T0 /= (self.r/delta)
        T2 /= (self.r/delta)
        cdef double predv = np.imag(1 + T0 * (self.r/delta)**3) * 1000*4/(MU0 * self.omega * self.r * self.r)
        cdef double predh = np.imag(1 + T2 * (self.r/delta)**2) * 1000*4/(MU0 * self.omega * self.r * self.r)

        return predv, predh


    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def _r0(self, double[:] lmbda, double[:] sigma):

        cdef int j, jj, k
        cdef int lsize = lmbda.size
        cdef int nlayer = self.g.nlayer

        cdef double[::1] d = self.g.d

        _Y = np.empty((nlayer-1), dtype=complex)
        _r0 = np.empty((lmbda.size,), dtype=complex)
        cdef double complex[::1] Y = _Y
        cdef double complex[::1] r0 = _r0

        _imo = 1j * self.g.mu * self.omega
        cdef double complex[::1] imo = _imo
        _ismo = np.empty((imo.size,), dtype=complex)
        cdef double complex[::1] ismo = _ismo
        for j in range(imo.size):
            ismo[j] = sigma[j] * imo[j]

        _u = np.empty((imo.size,), dtype=complex)
        _N = np.empty((imo.size,), dtype=complex)
        cdef double complex[::1] u = _u
        cdef double complex[::1] N = _N
        cdef double complex N0
        cdef double complex tanhud

        cdef double sigma0 = self.sigma0
        cdef double omega = self.omega

        for j in range(lsize):
            for jj in range(ismo.size):
                u[jj] = np.sqrt(lmbda[j]**2 + ismo[jj])
                N[jj] = u[jj] / imo[jj]

            tanhud = tanh(d[nlayer-2] * u[nlayer-2])
            Y[nlayer-2] = N[nlayer-2] * (N[nlayer-1] + N[nlayer-2] * tanhud) / \
                                        (N[nlayer-2] + N[nlayer-1] * tanhud)
            for k in range(nlayer-3, -1, -1):
                tanhud = tanh(d[k] * u[k])
                Y[k] = N[k] * (Y[k+1] + N[k] * tanhud) / (N[k] + Y[k+1] * tanhud)

            N0 = np.sqrt(lmbda[j]**2 + 1j * sigma0 * MU0 * omega) / (1j * MU0 * omega)
            r0[j] = (N0-Y[0]) / (N0+Y[0])
        return r0
