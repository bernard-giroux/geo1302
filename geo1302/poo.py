#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code utilitaire pour le cours GEO1302 - Partie Programmation orientée objet

Created on Wed Jan 18 14:35:23 2017

@author: giroux
"""
from abc import ABCMeta, abstractmethod

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class MaillageTriangulaire:
    def __init__(self, n, t):
        """Constructeur

        INPUT
            n: coordonnées X Z des noeuds (array n_noeuds x 2)
            t: indices des noeuds formant les triangles (array n_triangles x 3)
        """
        self.noeuds = n
        self.triangles = t

    def getLimites(self):
        """Retourne les min et max des noeuds
        """
        return np.hstack((self.noeuds.min(axis=0),
                          self.noeuds.max(axis=0)))

    @property
    def noeuds(self):
        return self._noeuds
    @noeuds.setter
    def noeuds(self, val):
        tmp = np.array(val, dtype=np.float64)
        if tmp.ndim != 2:
            raise ValueError('2D array needed')
        if tmp.shape[0] < 3:
            raise ValueError('3 nodes or more needed')
        if tmp.shape[1] != 2:
            raise ValueError('nodes: x z needed')

        self._noeuds = tmp

    @property
    def triangles(self):
        return self._triangles
    @triangles.setter
    def triangles(self, val):
        tmp = np.array(val, dtype=np.int)
        if tmp.ndim != 2:
            raise ValueError('2D array needed')
        if tmp.shape[1] != 3:
            raise ValueError('3 indices needed')

        self._triangles = tmp

    def __str__(self):
        return self.__class__.__name__+' ('+str(self.noeuds.shape[0])+' noeuds, '\
                                          +str(self.triangles.shape[0])+' triangles)'




class Donnees(metaclass=ABCMeta):
    def __init__(self, dtype):
        self.type = dtype

    @abstractmethod
    def getValeurs(self):
        pass
    @abstractmethod
    def getLocalisation(self):
        pass
    @abstractmethod
    def plot(self, fig=None):
        pass

class DonneesTerrain(Donnees):
    def __init__(self, dtype, data, loc):
        Donnees.__init__(self, dtype)
        self.val = np.array(data)
        self.locali = np.array(loc)

    def getValeurs(self):
        return self.val

    def getLocalisation(self):
        return self.locali

    def plot(self, fig=None):
        if fig == None:
            fig = plt.figure()

        norm = mpl.colors.Normalize(vmin=np.min(self.val), vmax=np.max(self.val))
        cmap = cm.jet
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        for n in np.arange(len(self.val)):
            plt.plot(self.locali[n,0], self.locali[n,1], 'o',
                     markerfacecolor=m.to_rgba(self.val[n]),
                     markeredgecolor=m.to_rgba(self.val[n]))
        return fig

class DonneesLabo(Donnees):
    def __init__(self, dtype, val, loc, no):
        Donnees.__init__(self, dtype)
        self.val = val
        self.loc = loc
        self.no_ech = no

    def getValeurs(self):
        return self.val

    def getLocalisation(self):
        return self.loc

    def plot(self, fig=None):
        if fig == None:
            fig = plt.figure()

        plt.plot(self.no_ech, self.val, 'sk')
        return fig

if __name__ == '__main__':

    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    tri = np.array([[0, 1, 2]])
    mt = MaillageTriangulaire(nodes, tri)

#    tri = np.array([[0, 1, 2, 3]])
#    mt2 = MaillageTriangulaire(nodes, tri);

    mt2 = mt
    print(mt.triangles)
    mt2.triangles = np.array([[5, 6, 7]])
    print(mt.triangles)

    import copy

    mt2 = copy.copy(mt)
    mt2.triangles = np.array([[3, 4, 5]])
    print(mt.triangles)
    print(mt2.triangles)

    class ID:
        def __init__(self, v = 0):
            self.value = v

    mt.id = ID()
    mt2 = copy.copy(mt)
    print(mt.id.value)
    print(mt2.id.value)
    mt2.id.value = 2
    print(mt.id.value)
    print(mt2.id.value)



    mt2 = copy.deepcopy(mt)
    print(mt.id.value)
    print(mt2.id.value)
    mt2.id.value = 4
    print(mt.id.value)
    print(mt2.id.value)


    print(mt)

    dt = DonneesTerrain('gravi', [4.0, 6.4, 3.1], [[3.0, 4.5],[7.6, 3.2],[6.5, 4.2]])
    dt.plot()

    dl = DonneesLabo('densité', 2.67, 'Site no 2', 1)
    f = dl.plot()
    dl2 = DonneesLabo('densité', 2.62, 'Site no 2', 2)
    dl2.plot(f)
