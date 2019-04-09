#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 22:33:51 2019

@author: giroux
"""
import numpy as np


from em import EM38

def get_l_rough(n, deg):
    if deg == 0:
        return np.eye(n)

    df = np.hstack((np.array([-1, 1]), np.zeros((deg-1,))))
    for i in range(1, deg):
        # take the difference of the lower order derivative and itself shifted left
        # to get a derivative one order higher
        df = np.hstack((np.zeros((1,)), df[:deg])) - np.hstack((df[:deg], np.zeros((1,))))

    dn  = n - deg
    D = np.zeros((dn, n))
    for d in range(deg+1):
        for i in range(dn):
            j = i+d
            D[i, j] = df[d]
    return D


def occam(G, jac, D, d, m0, delta):

    m = m0.copy()
    oldm = np.zeros(m.shape)
    it = 0

    DTD = D.T.dot(D)

    Gm = G(m)
    mchi2 = np.linalg.norm(Gm - d)**2
    alphas = np.logspace(-20, 0, 100)
    chis = np.empty((alphas.size,))

    while (np.linalg.norm(oldm-m)/np.linalg.norm(m) > 5.0e-3) or (mchi2 > delta*delta*1.01):
        it += 1
        if it > 30:
            return m

        oldm = m.copy()

        J = jac(m)

        dhat =

        for i in range(alphas.size):
            M =
            if np.linalg.cond(M) < 1.e15:

                m = np.linalg.solve(M, J.T.dot(dhat))

                chis[i] = np.linalg.norm(G(m) - d)**2
            else:
                chis[i] = np.inf

        y = chis.min()
        if y > delta*delta:
            i = np.argmin(chis)
            alpha = alphas[i]
        else:
            i = np.where(chis<=delta*delta)
            alpha = alphas[np.max(i)]

        m = np.linalg.solve(
        Gm = G(m)
        mchi2 = np.linalg.norm(Gm - d)**2

    return m


if __name__ == '__main__':

    heights = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.50])
    mtrue = np.array([100, 95, 90, 95, 100, 130, 160, 200, 250, 310, 360])/1000.0


    em38 = EM38(heights)

    datanf = em38.G(mtrue)
    data = datanf + 0.1*np.random.randn(datanf.size)

    m0 = 200.0/1000.0 + np.zeros(mtrue.shape)

    D = get_l_rough(mtrue.size, 2)
    delta = 0.1 * np.sqrt(18)
    moccam = occam(em38.G, em38.jac, D, data, m0, delta)
