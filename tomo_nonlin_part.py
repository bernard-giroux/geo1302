#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 12:58:03 2018

@author: giroux
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.sparse.linalg as spl

import ttcrpy.rgrid as rg

plt.style.use('default')


# %% fonction pour l'affichage des rais

def plot_rais(rais, grid=None, ax=None):
    """

    Affiche les rais

    Parameters
    ----------
    rais : :obj:`list` of :obj:`np.ndarray`
        Coordonnées des rais.
    grid : :obj:`Grid2d`, optional
        Grille de modélisation pour l'affichage du modèle de vitesse. The
        default is None.
    ax : :obj:`Axes`, optional
        Axes où tracer les rais, llLes axes sont créés automatiquement si
        None. The default is None.

    Returns
    -------
    ax : :obj:`Axes`
        Axes utilisés pour tracer les rais.

    """
    if ax is None:
        plt.figure(figsize=(6,8))
        ax = plt.gca()

    if grid is not None:
        nx, nz = grid.shape
        x = grid.x
        z = grid.z
        slowness = grid.get_slowness()

        im = ax.pcolor(x, z, 1/(slowness.T), cmap='plasma')
        cb = plt.colorbar(im)
        cb.set_label('Vitesse (m/s)')

    for ir in range(len(rais)):
        ax.plot(rais[ir][:, 0], rais[ir][:, 1], 'k', linewidth=0.5)

    plt.axis('scaled')
    ax.invert_yaxis()
    ax.set_ylabel('Profondeur (m)')
    ax.set_xlabel('Distance (m)')
    return ax


# %% Construction de la grille pour le tracé de rais

dx = 5.0
dz = 5.0
xmin = 0.0
zmin = 0.0
xmax = 150.0
zmax = 300.0
nx = int((xmax-xmin)/dx+0.001)
nz = int((zmax-xmin)/dz+0.001)
xg = np.linspace(xmin, xmax, nx+1)
zg = np.linspace(zmin, zmax, nz+1)

g = rg.Grid2d(xg, zg, method='SPM', nsnx=10, nsnz=10)

# %% Lecture des données

data = np.loadtxt('model1_tt.dat')
Tx = data[:, :2].copy()
Rx = data[:, 2:4].copy()
dobs = data[:, -1]
t0 = np.zeros((Tx.shape[0], ))


# %%
# modèle initial, on prend la lenteur apparente moyenne
Ldroit = rg.Grid2d.data_kernel_straight_rays(Tx, Rx, xg, zg)  # matrice L pour des rais droits
l_rai =                 # longueur des rais
lent_app = dobs/l_rai                     # lenteurs apparentes
m0 = np.mean(lent_app) + np.zeros((nx*nz,))


# %%

# matrices de lissage
Dx, Dz = g.compute_K(2)

D = Dx + Dz
DTD = D.T.dot(D)


# %%
alpha_val = [10.0, 15.0, 25.0, 40.0, 80.0, 160.0, 300.0, 500.0, 1000.0, 10000.]
maxit = 4
nfig = 1

models = []
residuals = []
sol_norm = []
for alpha in alpha_val:
    E = []
    for it in range(maxit):
        if it == 0:
            J = Ldroit
            m = m0.copy()
            d = Ldroit.dot(m0)
        else:
            d, L = g.raytrace(Tx, Rx, slowness=m, compute_L=True)
            J = L

        E.append(np.linalg.norm(d-dobs))
        print('Itération {0:d}, E = {1:5.4e}'.format(it+1, E[it]))

        A =
        b =
        x = spl.lsqr(A, b)
        dm = x[0]
        m += dm

        if it == maxit-1:
            d, L = g.raytrace(Tx, Rx, slowness=m, compute_L=True)
            E.append(np.linalg.norm(d-dobs))

            plt.figure(nfig, figsize=(9, 5))
            plt.subplot(121)
            plt.pcolor(xg, zg, 1/(m.reshape(nx, nz).T), cmap='plasma')
            plt.clim(2100.0, 2800.0)
            plt.axis('scaled')
            plt.gca().invert_yaxis()
            plt.ylabel('Profondeur (m)', fontsize=12)
            plt.xlabel('Distance (m)', fontsize=12)
            cb = plt.colorbar()

            cb.set_label('Vitesse (m/s)', fontsize=12)
            plt.subplot(122)
            plt.plot(np.arange(maxit+1), E)
            plt.ylabel('E', fontsize=12)
            plt.xlabel('Itération', fontsize=12)
            plt.suptitle('α = {0:3.1f}'.format(alpha), fontsize=14)
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            plt.savefig('fig/tomo_m1_a{0:d}.pdf'.format(int(alpha)), bbox_inches='tight')
            plt.show(block=False)
            plt.pause(1)

            models.append(m)
            residuals.append(E[-1])
            sol_norm.append(np.linalg.norm(D.dot(m)))

    nfig += 1

# %%

plt.style.use('seaborn')

plt.figure(nfig, figsize=(6, 5))
plt.plot(residuals, sol_norm, 'k-', zorder=1)
plt.scatter(residuals, sol_norm, c=alpha_val, cmap='plasma', zorder=2,
            norm=colors.LogNorm(vmin=alpha_val[0], vmax=alpha_val[-1]), s=50)
plt.xlabel('$\|\|G(\mathbf{m})-\mathbf{d}_{obs}\|\|_2$', fontsize=16)
plt.ylabel('$\|\|\mathbf{Dm}\|\|_2$', fontsize=16)
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
cb = plt.colorbar()
cb.ax.set_ylabel('$\\alpha$', fontsize=16)
plt.tight_layout()
plt.savefig('fig/tomo_m1_norm.pdf'.format(alpha))
plt.show()
