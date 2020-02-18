#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 22:25:43 2017

@author: giroux
"""
import math
import scipy.sparse as sp
from scipy.interpolate import interpn
from scipy.sparse.linalg import bicgstab
import numpy as np
import h5py

from geo1302.utils import nargout


class GrilleVF:
    """
    Classe pour gérer des grilles régulières
    """
    def __init__(self, x, y, z):
        """
        Input
            x: coordonnées des noeuds selon x
            y: coordonnées des noeuds selon y
            z: coordonnées des noeuds selon z
        """
        self.x = x
        self.y = y
        self.z = z
        self.__update_n()

    def ind(self, i, j, k):
        """
        Retourne les indices de prismes

        Input
            i: indice(s) du voxel selon x
            j: indice(s) du voxel selon y
            k: indice(s) du voxel selon z
        """
        if np.size(i)>1:
            i = np.array(i)
            i = i.flatten()
        if np.size(j)>1:
            j = np.array(j)
            j = j.flatten()
        if np.size(k)>1:
            k = np.array(k)
            k = k.flatten()
        if np.any(i<0) or
            raise IndexError('Index out of bound')
        # on utlise np.size car i peut être un scalaire
        if np.size(i)>1 or np.size(j)>1 or np.size(k)>1:
            ii = np.kron(
            jj = np.kron(
            kk = np.kron(
            return np.sort((kk * self.ny + jj) * self.nx + ii)
        else:
            return (k * self.ny + j) * self.nx + i

    def __update_n(self):
        if '_x' in self.__dict__ and '_y' in self.__dict__ and '_z' in self.__dict__:
            self.nx = len(self._x)-1             # nombre de prismes selon x
            self.ny = len(self._y)-1             # nombre de prismes selon y
            self.nz = len(self._z)-1             # nombre de prismes selon z
            self.nc = self.nx * self.ny * self.nz    # nombre de prismes
            self.nfx =
            self.nfy =
            self.nfz =
            self.nf =


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
            if len(tmp)<2:
                raise ValueError('2 nodes or more needed')
        except:
            raise

        self._x = tmp
        self.hx = np.diff(tmp)
        self.xc = (tmp[1:]+tmp[:-1])/2
        self.dx = np.diff(self.xc)
        self.__update_n()

    @property
    def y(self):
        "Coordonnées des noeuds selon y"
        return self._y

    @y.setter
    def y(self, val):
        try:
            tmp = np.array(val, dtype=np.float64)
            if tmp.ndim != 1:
                raise ValueError('1D array needed')
            if len(tmp)<2:
                raise ValueError('2 nodes or more needed')
        except:
            raise

        self._y = tmp
        self.hy = np.diff(tmp)
        self.yc = (tmp[1:]+tmp[:-1])/2
        self.dy = np.diff(self.yc)
        self.__update_n()

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
            if len(tmp)<2:
                raise ValueError('2 nodes or more needed')
        except:
            raise

        self._z = tmp
        self.hz = np.diff(tmp)
        self.zc = (tmp[1:]+tmp[:-1])/2
        self.dz = np.diff(self.zc)
        self.__update_n()


    def magmod(self, chi, B0, xo, spheq, chtot, tol=1.e-8, maxit=1000,
               solver=bicgstab, precon=False):
        '''
        Calcul la réponse magnétique

        Paramètres
        ----------
        chi    : susceptibilité des voxels de la grille (nc x 1)
        B0     : champ ambiant (3 x 1)
        xo     : pts d'observation (npts x 3)
        spheq  : sphère équivalente pour les cond. limites (1: oui, 0:non)
        chtot  : calcul le champ total (1) ou l'anomalie (0)
        tol    : tolérance du solveur
        maxit  : Nbre max d'itération du solveur
        solver : solveur itératif du module scipy.sparse.linalg
        precon : si True, utilise un préconditionneur (factorisation LU incomplète)

        Retourne
        --------
        Bx, By, Bz  : arrays numpy
            valeurs du champ interpolé au points xo

        '''
        B0 = np.array(B0)
        mu0 = 4*math.pi*1.e-7
        mu = mu0 * (1+chi)

        D = self.fabrique_D()
        M = self.fabrique_M(mu)
        G = self.fabrique_G()
        f,g = self.fabrique_cf(B0,chi)

        if spheq == 0:
            g[:] = 0.0

        if chtot == 1:
            q = f + g
            A = D.dot(M.dot(G))
            if precon:
                Mpre = sp.linalg.spilu(A.tocsc())
                Mpre = sp.linalg.LinearOperator(A.shape, Mpre.solve)
            else:
                Mpre = None
            phi, info = solver(A, q, tol=tol, maxiter=maxit, M=Mpre)
            B = M.dot(G.dot(phi))
        else:
            BB0 = np.hstack((B0[0]+np.zeros((self.nfx,)),
                             B0[1]+np.zeros((self.nfy,)),
                             B0[2]+np.zeros((self.nfz,))))

            Mtmp =
            q =
            A =
            if precon:
                Mpre = sp.linalg.spilu(A.tocsc())
                Mpre = sp.linalg.LinearOperator(A.shape, Mpre.solve)
            else:
                Mpre = None
            phi_s, info = solver(A, q, tol=tol, maxiter=maxit, M=Mpre)
            B =

        if info > 0:
            print('{0:}: convergence not achieved, stopped after {1:d} iterations for tol = {2:g}'.format(solver.__name__, info, tol))
        elif info < 0:
            print('{0:s}: illegal input or breakdown'.format(solver.__name__))

        Bx = B[:self.nfx]
        By = B[self.nfx:(self.nfx+self.nfy)]
        Bz = B[(self.nfx+self.nfy):]

        Bx = Bx.reshape(self.nx-1,self.ny,self.nz, order='F')
        By = By.reshape(self.nx,self.ny-1,self.nz, order='F')
        Bz = Bz.reshape(self.nx,self.ny,self.nz-1, order='F')

        Bx = interpn((self.x[1:-1], self.yc, self.zc), Bx, xo)
        By = interpn((self.xc, self.y[1:-1], self.zc), By, xo)
        Bz = interpn((self.xc, self.yc, self.z[1:-1]), Bz, xo)

        return Bx, By, Bz

    def sol_an_dipole(self, chi, B0, xo):
        """
        Solution analytique pour un dipole équivalent
        """

        # réplique des largeurs des voxels pour l'ensemble des voxels
        hx = np.kron(np.ones((self.ny*self.nz,)), self.hx)
        hy = np.kron(np.kron(np.ones((self.nz,)), self.hy), np.ones((self.nx,)))
        hz = np.kron(self.hz, np.ones(self.nx*self.ny,))

        # réplique des coord des centres des voxels pour l'ensemble des voxels
        xc = np.kron(np.ones((self.ny*self.nz,)), self.xc)
        yc = np.kron(np.kron(np.ones((self.nz,)), self.yc), np.ones((self.nx,)))
        zc = np.kron(self.zc, np.ones((self.nx*self.ny,)))

        ind = np.logical_not( chi==0.0 )  # voxels magnétisables
        v = hx[ind] * hy[ind] * hz[ind]
        V = np.sum(v)
        xi = np.sum( chi[ind] * v )/V  # éq. 4-1b de Lelièvre (2003)
        m = B0 / (4*math.pi*1.e-7) * V*xi/(1.+xi/3.) # éq. 4-5 de Lelièvre (2003)
        mm = np.linalg.norm(m)  # magnitude
        m = m.reshape(-1,1)/mm                # vecteur unitaire

        # "centre de susceptibilité"
        xc = np.sum(chi[ind] * xc[ind]) / np.sum(chi[ind])
        yc = np.sum(chi[ind] * yc[ind]) / np.sum(chi[ind])
        zc = np.sum(chi[ind] * zc[ind]) / np.sum(chi[ind])

        r = np.vstack((xo[:,0]-xc, xo[:,1]-yc, xo[:,2]-zc))
        rm = np.sqrt( np.sum(r*r, axis=0))
        r = r/np.kron(rm, np.ones((3,1)))     # vecteur unitaire

        mtmp = np.kron(m, np.ones((rm.size,)))
        mtmp = np.kron(3*np.sum(mtmp*r,axis=0),np.ones((3,1)))*r - mtmp
        B = np.kron(1.e-7*mm/(rm*rm*rm), np.ones((3,1))) * mtmp

        return B[0,:], B[1,:], B[2,:]

    def fabrique_D(self):

        # Dx

        M = self.nx
        N = self.nx-1
        i =
        j =
        nval = i.size
        ii = np.empty((self.ny*self.nz*nval,), dtype=np.int64)
        jj = np.empty((self.ny*self.nz*nval,), dtype=np.int64)
        for n in np.arange(self.ny*self.nz):
            ii[] =
            jj[] =
        s = np.tile(np.hstack((1.0/self.hx[:-1], -1.0/self.hx[1:])), (self.ny*self.nz,))
        Dx = sp.coo_matrix((s,(ii,jj)))

        # Dy

        M = self.nx*self.ny
        N = self.nx*(self.ny-1)
        i =
        j =
        nval = i.size
        ii = np.empty((self.nz*nval,), dtype=np.int64)
        jj = np.empty((self.nz*nval,), dtype=np.int64)
        for n in np.arange(self.nz):
            ii[] =
            jj[] =
        s = np.hstack((np.kron(1.0/self.hy[:-1], np.ones((self.nx,))),
                       np.kron(-1.0/self.hy[1:], np.ones((self.nx,)))))
        s = np.tile(s, (self.nz,))
        Dy = sp.coo_matrix((s, (ii,jj)))

        # Dz

        N = (self.nx*self.ny*(self.nz-1))
        i =
        j =
        s = np.hstack((np.kron(1./self.hz[:-1], np.ones((self.nx*self.ny,))),
                       np.kron(-1./self.hz[1:], np.ones((self.nx*self.ny,))) ))
        Dz = sp.coo_matrix((s, (i,j)))

        # assemblage

        return sp.hstack((Dx, Dy, Dz)).tocsr()

    def fabrique_q(self, B0):
        """
        Création du vecteur des conditions aux frontières

        Input
        B0: champ ambiant assigné aux frontières (vecteur [Bx, By, Bz])
        """
        q = np.zeros((self.nc,))

        # Face correspondant à x_min
        ind = self.ind(0, np.arange(self.ny), np.arange(self.nz))
        q[ind] = B0[0]/self.hx[0]

        # Face correspondant à x_max
        ind =
        q[ind] = -B0[0]/self.hx[-1]


        # Face correspondant à y_min
        ind =
        q[ind] = q[ind] + B0[1]/self.hy[0]

        # Face correspondant à y_max
        ind =
        q[ind] = q[ind] - B0[1]/self.hy[-1]


        # Face correspondant à z_min
        ind =
        q[ind] = q[ind] + B0[2]/self.hz[0]

        # Face correspondant à z_max
        ind =
        q[ind] = q[ind] - B0[2]/self.hz[-1]

        return q

    def fabrique_M(self, mu):

        # Mx

        d =
        hi = self.hx[1:]
        him1 = self.hx[:-1]
        tmp = np.ones((self.nx,), dtype=bool)
        tmp[0] = 0
        ind1 = np.kron(np.ones((self.ny*self.nz,), dtype=bool), tmp)
        tmp[0] = 1
        tmp[-1] = 0
        ind2 =
        Mx = d / ( hi/mu[ind1] + him1/mu[ind2] )

        # My

        d =
        hi = np.kron( self.hy[1:],
        him1 = np.kron( self.hy[:-1],
        tmp = np.ones((self.ny,), dtype=bool)
        tmp[0] = 0
        ind1 =
        tmp[0] = 1
        tmp[-1] = 0
        ind2 =
        My = d / ( hi/mu[ind1] + him1/mu[ind2] )

        # Mz

        d =
        hi = np.kron(self.hz[1:], np.ones((self.nx*self.ny,)))
        him1 = np.kron(self.hz[:-1], np.ones((self.nx*self.ny,)))
        ind1 =
        ind2 = np.hstack((np.ones((self.nx*self.ny*(self.nz-1),), dtype=bool),
                          np.zeros((self.nx*self.ny,), dtype=bool) ))
        Mz = d / ( hi/mu[ind1] + him1/mu[ind2] )

        # assemblage

        return sp.coo_matrix((np.hstack((Mx, My, Mz)),
                                        (np.arange(self.nf),
                                         np.arange(self.nf)))).tocsr()

    def fabrique_G(self):

        # Gx

        M = self.nx - 1
        N = self.nx
        i =
        j =
        nval = i.size
        ii = np.empty((self.ny*self.nz*nval,), dtype=np.int64)
        jj = np.empty((self.ny*self.nz*nval,), dtype=np.int64)
        for n in np.arange(self.ny*self.nz):
            ii[] =
            jj[] =
        s = np.tile(np.hstack((-1.0/self.dx, 1.0/self.dx)), (self.ny*self.nz,))
        Gx = sp.coo_matrix((s,(ii,jj)))

        # Gy

        M = self.nx*(self.ny-1)
        N = self.nx*self.ny
        i =
        j =
        nval = i.size
        ii = np.empty((self.nz*nval,), dtype=np.int64)
        jj = np.empty((self.nz*nval,), dtype=np.int64)
        for n in np.arange(self.nz):
            ii[] =
            jj[] =
        s = np.hstack((np.kron(-1.0/self.dy, np.ones((self.nx,))),
                       np.kron( 1.0/self.dy, np.ones((self.nx,)))))
        s = np.tile(s, (self.nz,))
        Gy = sp.coo_matrix((s,(ii,jj)))

        # Gz

        i =
        j =
        s = np.hstack((np.kron(-1.0/self.dz, np.ones((self.nx*self.ny,))),
                       np.kron( 1.0/self.dz, np.ones((self.nx*self.ny,))) ))
        Gz = sp.coo_matrix((s,(i,j)))

        # assemblage
        return sp.vstack((Gx, Gy, Gz)).tocsr()

    def fabrique_cf(self, B0, chi=[]):
        nout = nargout()

        if not ( nout==1 or nout==2 ):
            raise ValueError()

        if len(chi)==0 and nout > 1:
            raise ValueError('Erreur: donnez chi en argument pour obtenir f et g')

        f = np.zeros((self.nc))

        if len(chi)>0:
            g = np.zeros((self.nc))

            # réplique des largeurs des voxels pour l'ensemble des voxels
            hx = np.kron(np.ones((self.ny*self.nz,)),self.hx)
            hy = np.kron(np.kron(np.ones((self.nz,)),self.hy),np.ones((self.nx,)))
            hz = np.kron(self.hz,np.ones((self.nx*self.ny,)))

            # réplique des coord des centres des voxels pour l'ensemble des voxels
            xc = np.kron(np.ones((self.ny*self.nz,)),self.xc)
            yc = np.kron(np.kron(np.ones((self.nz,)),self.yc),np.ones((self.nx,)))
            zc = np.kron(self.zc,np.ones((self.nx*self.ny,)))

            ind = np.logical_not( chi==0.0 ) # voxels magnétisables
            v = hx[ind] * hy[ind] * hz[ind]
            V = np.sum(v)
            xi = np.sum( chi[ind] * v )/V  # éq. 4-1b de Lelièvre (2003)
            m = B0 / (4*math.pi*1.e-7) * V*xi/(1.+xi/3.) # éq. 4-5 de Lelièvre (2003)
            mm = np.linalg.norm(m)  # magnitude
            m = m.reshape(-1,1)/mm                # vecteur unitaire

            # "centre de susceptibilité"
            xc = np.sum(chi[ind] * xc[ind]) / np.sum(chi[ind])
            yc = np.sum(chi[ind] * yc[ind]) / np.sum(chi[ind])
            zc = np.sum(chi[ind] * zc[ind]) / np.sum(chi[ind])

        # -------------------------------------------------------------
        # X
        # -------------------------------------------------------------

        # Face correspondant à x_min
        ind = self.ind(0, np.arange(self.ny), np.arange(self.nz))
        f[ind] = B0[0] / self.hx[0]
        if len(chi)>0:
            # coord des pts d'observation
            xp = self.x[0] + np.zeros((self.ny*self.nz,))
            yp =
            zp = np.kron(self.zc, np.ones((self.ny,)))

            r = np.vstack((xp-xc, yp-yc, zp-zc))
            rm = np.sqrt( np.sum(r*r, axis=0))
            ind2 = rm>0
            r = r/np.kron(rm, np.ones((3,1)))   # vecteur unitaire

            mtmp = np.kron(m, np.ones((rm.size,)))
            mtmp = np.kron(3*np.sum(mtmp*r,axis=0),np.ones((3,1)))*r - mtmp
            B = np.kron(1.e-7*mm/(rm*rm*rm), np.ones((3,1))) * mtmp
            g[ind[ind2]] = B[0,ind2] / self.hx[0]

        # Face correspondant à x_max
        ind =
        f[ind] = -B0[0] / self.hx[-1]
        if len(chi)>0:
            # coord des pts d'observation
            xp =
            yp =
            zp =

            r = np.vstack((xp-xc, yp-yc, zp-zc))
            rm = np.sqrt( np.sum(r*r, axis=0))
            ind2 = rm>0
            r = r/np.kron(rm, np.ones((3,1)))

            mtmp = np.kron(m, np.ones((rm.size,)))
            mtmp = np.kron(3*np.sum(mtmp*r,axis=0),np.ones((3,1)))*r - mtmp
            B = np.kron(1.e-7*mm/(rm*rm*rm), np.ones((3,1))) * mtmp
            g[ind[ind2]] = -B[0,ind2] / self.hx[-1]

        # -------------------------------------------------------------
        # Y
        # -------------------------------------------------------------

        # Face correspondant à y_min
        ind =
        f[ind] += B0[1] / self.hy[0]
        if len(chi)>0:
            # coord des pts d'observation
            xp =
            yp =
            zp =

            r = np.vstack((xp-xc, yp-yc, zp-zc))
            rm = np.sqrt( np.sum(r*r, axis=0))
            ind2 = rm>0
            r = r/np.kron(rm, np.ones((3,1)))

            mtmp = np.kron(m, np.ones((rm.size,)))
            mtmp = np.kron(3*np.sum(mtmp*r,axis=0),np.ones((3,1)))*r - mtmp
            B = np.kron(1.e-7*mm/(rm*rm*rm), np.ones((3,1))) * mtmp
            g[ind[ind2]] += B[1,ind2] / self.hy[0]

        # Face correspondant à y_max
        ind =
        f[ind] -= B0[1] / self.hy[-1]
        if len(chi)>0:
            # coord des pts d'observation
            xp =
            yp =
            zp =

            r = np.vstack((xp-xc, yp-yc, zp-zc))
            rm = np.sqrt( np.sum(r*r, axis=0))
            ind2 = rm>0
            r = r/np.kron(rm, np.ones((3,1)))

            mtmp = np.kron(m, np.ones((rm.size,)))
            mtmp = np.kron(3*np.sum(mtmp*r,axis=0),np.ones((3,1)))*r - mtmp
            B = np.kron(1.e-7*mm/(rm*rm*rm), np.ones((3,1))) * mtmp
            g[ind[ind2]] -= B[1,ind2] / self.hy[-1]

        # -------------------------------------------------------------
        # Z
        # -------------------------------------------------------------

        # Face correspondant à z_min
        ind =
        f[ind] += B0[2] / self.hz[0]
        if len(chi)>0:
            # coord des pts d'observation
            xp =
            yp =
            zp =

            r = np.vstack((xp-xc, yp-yc, zp-zc))
            rm = np.sqrt( np.sum(r*r, axis=0))
            ind2 = rm>0
            r = r/np.kron(rm, np.ones((3,1)))

            mtmp = np.kron(m, np.ones((rm.size,)))
            mtmp = np.kron(3*np.sum(mtmp*r,axis=0),np.ones((3,1)))*r - mtmp
            B = np.kron(1.e-7*mm/(rm*rm*rm), np.ones((3,1))) * mtmp
            g[ind[ind2]] += B[2,ind2] / self.hz[0]

        # Face correspondant à z_max
        ind =
        f[ind] -= B0[2] / self.hz[-1]
        if len(chi)>0:
            # coord des pts d'observation
            xp =
            yp =
            zp =

            r = np.vstack((xp-xc, yp-yc, zp-zc))
            rm = np.sqrt( np.sum(r*r, axis=0))
            ind2 = rm>0
            r = r/np.kron(rm, np.ones((3,1)))

            mtmp = np.kron(m, np.ones((rm.size,)))
            mtmp = np.kron(3*np.sum(mtmp*r,axis=0),np.ones((3,1)))*r - mtmp
            B = np.kron(1.e-7*mm/(rm*rm*rm), np.ones((3,1))) * mtmp
            g[ind[ind2]] -= B[2,ind2] / self.hz[-1]

        # Output

        if nout==1 and len(chi)==0:
            return f
        elif nout==2 and len(chi)>0:
            return f, g
        else: # nécessairement nout==1 && len(chi)>0
            return f+g

    def toXdmf(self, field, fieldname, filename):
        """
        Save a field in xdmf format (http://www.xdmf.org/index.php/Main_Page)

        INPUT
            field: data array of size equal to the number of cells in the grid
            fieldname: name to be assigned to the data (string)
            filename: name of xdmf file (string)
        """

        f = open(filename+'.xmf','w')

        f.write('<?xml version="1.0" ?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        f.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
        f.write(' <Domain>\n')
        f.write('   <Grid Name="Structured Grid" GridType="Uniform">\n')
        f.write('     <Topology TopologyType="3DRECTMesh" NumberOfElements="'+
        repr(self.nz+1)+' '+repr(self.ny+1)+' '+repr(self.nx+1)+'"/>\n')
        f.write('     <Geometry GeometryType="VXVYVZ">\n')
        f.write('       <DataItem Dimensions="'+repr(self.nz)+'" NumberType="Float" Precision="4" Format="XML">\n')
        for z in self.z:
            f.write('          '+repr(z)+'\n')
        f.write('       </DataItem>\n')
        f.write('       <DataItem Dimensions="'+repr(self.ny)+'" NumberType="Float" Precision="4" Format="XML">\n')
        for y in self.y:
            f.write('          '+repr(y)+'\n')
        f.write('       </DataItem>\n')
        f.write('       <DataItem Dimensions="'+repr(self.nx)+'" NumberType="Float" Precision="4" Format="XML">\n')
        for x in self.x:
            f.write('          '+repr(x)+'\n')
        f.write('       </DataItem>\n')
        f.write('     </Geometry>\n')
        f.write('     <Attribute Name="'+fieldname+'" AttributeType="Scalar" Center="Cell">\n')
        f.write('       <DataItem Dimensions="'+repr(self.nz)+' '+repr(self.ny)+' '+repr(self.nx)+
        '" NumberType="Float" Precision="4" Format="HDF">'+filename+'.h5:/'+fieldname+'</DataItem>\n')
        f.write('     </Attribute>\n')
        f.write('   </Grid>\n')
        f.write(' </Domain>\n')
        f.write('</Xdmf>\n')

        f.close()

        h5f = h5py.File(filename+'.h5', 'w')
        h5f.create_dataset(fieldname, data=field.reshape((self.nz,self.ny,self.nx)).astype(np.float32))
        h5f.close()

    def toVTK(self, field, fieldname, filename):
        """
        Save a field in VTK format

        INPUT
            field: data array of size equal to the number of cells in the grid
            fieldname: name to be assigned to the data (string)
            filename: name of vtk file without extension (string)
        """
        import vtk

        xCoords = vtk.vtkFloatArray()
        for i in self.x:
            xCoords.InsertNextValue(i)

        yCoords = vtk.vtkFloatArray()
        for i in self.y:
            yCoords.InsertNextValue(i)

        zCoords = vtk.vtkFloatArray()
        for i in self.z:
            zCoords.InsertNextValue(i)

        rgrid = vtk.vtkRectilinearGrid()
        rgrid.SetDimensions(len(self.x), len(self.y), len(self.z))
        rgrid.SetXCoordinates(xCoords)
        rgrid.SetYCoordinates(yCoords)
        rgrid.SetZCoordinates(zCoords)

        data = vtk.vtkDoubleArray()
        for i in field:
            data.InsertNextValue(i)
        data.SetName(fieldname)
        rgrid.GetCellData().AddArray(data)

        writer = vtk.vtkXMLRectilinearGridWriter()
        writer.SetInputData(rgrid)
        writer.SetFileName(filename+'.vtr')
        writer.Write()

    def toVTKx(self, field, fieldname, filename):
        """
        Save a x component field in VTK format

        INPUT
            field: data array of size equal to the number of cells in the grid
                   for the x component (nx-1) * ny * nz
            fieldname: name to be assigned to the data (string)
            filename: name of vtk file without extension (string)
        """
        import vtk

        xCoords = vtk.vtkFloatArray()
        for i in self.xc:
            xCoords.InsertNextValue(i)

        yCoords = vtk.vtkFloatArray()
        for i in self.y:
            yCoords.InsertNextValue(i)

        zCoords = vtk.vtkFloatArray()
        for i in self.z:
            zCoords.InsertNextValue(i)

        rgrid = vtk.vtkRectilinearGrid()
        rgrid.SetDimensions(len(self.xc), len(self.y), len(self.z))
        rgrid.SetXCoordinates(xCoords)
        rgrid.SetYCoordinates(yCoords)
        rgrid.SetZCoordinates(zCoords)

        data = vtk.vtkDoubleArray()
        for i in field:
            data.InsertNextValue(i)
        data.SetName(fieldname)
        rgrid.GetCellData().AddArray(data)

        writer = vtk.vtkXMLRectilinearGridWriter()
        writer.SetInputData(rgrid)
        writer.SetFileName(filename+'.vtr')
        writer.Write()

    def toVTKy(self, field, fieldname, filename):
        """
        Save a y component field in VTK format

        INPUT
            field: data array of size equal to the number of cells in the grid
                   for the y component nx * (ny-1) * nz
            fieldname: name to be assigned to the data (string)
            filename: name of vtk file without extension (string)
        """
        import vtk

        xCoords = vtk.vtkFloatArray()
        for i in self.x:
            xCoords.InsertNextValue(i)

        yCoords = vtk.vtkFloatArray()
        for i in self.yc:
            yCoords.InsertNextValue(i)

        zCoords = vtk.vtkFloatArray()
        for i in self.z:
            zCoords.InsertNextValue(i)

        rgrid = vtk.vtkRectilinearGrid()
        rgrid.SetDimensions(len(self.x), len(self.yc), len(self.z))
        rgrid.SetXCoordinates(xCoords)
        rgrid.SetYCoordinates(yCoords)
        rgrid.SetZCoordinates(zCoords)

        data = vtk.vtkDoubleArray()
        for i in field:
            data.InsertNextValue(i)
        data.SetName(fieldname)
        rgrid.GetCellData().AddArray(data)

        writer = vtk.vtkXMLRectilinearGridWriter()
        writer.SetInputData(rgrid)
        writer.SetFileName(filename+'.vtr')
        writer.Write()

    def toVTKz(self, field, fieldname, filename):
        """
        Save a z component field in VTK format

        INPUT
            field: data array of size equal to the number of cells in the grid
                   for the z component nx * ny * (nz-1)
            fieldname: name to be assigned to the data (string)
            filename: name of vtk file without extension (string)
        """
        import vtk

        xCoords = vtk.vtkFloatArray()
        for i in self.x:
            xCoords.InsertNextValue(i)

        yCoords = vtk.vtkFloatArray()
        for i in self.y:
            yCoords.InsertNextValue(i)

        zCoords = vtk.vtkFloatArray()
        for i in self.zc:
            zCoords.InsertNextValue(i)

        rgrid = vtk.vtkRectilinearGrid()
        rgrid.SetDimensions(len(self.x), len(self.y), len(self.zc))
        rgrid.SetXCoordinates(xCoords)
        rgrid.SetYCoordinates(yCoords)
        rgrid.SetZCoordinates(zCoords)

        data = vtk.vtkDoubleArray()
        for i in field:
            data.InsertNextValue(i)
        data.SetName(fieldname)
        rgrid.GetCellData().AddArray(data)

        writer = vtk.vtkXMLRectilinearGridWriter()
        writer.SetInputData(rgrid)
        writer.SetFileName(filename+'.vtr')
        writer.Write()


if __name__ == '__main__':

    x = [1, 2, 3, 3.5]
    y = [1, 2, 3, 4, 5]
    z = np.arange(6)
    gvf = GrilleVF(x, y, z)

    print(gvf.nc, gvf.dx)

    ind = gvf.ind([1,2],2,[3,0])
    print(ind)

    B0 = np.array([1., 2., 3.])
    chi = np.zeros((gvf.nc,))
    chi[gvf.ind(2,2,3)] = 1.0
    mu0 = 4 * math.pi * 1.e-7;
    mu = mu0 * (1.+chi)
    D = gvf.fabrique_D()
    M = gvf.fabrique_M(mu)
    G = gvf.fabrique_G()
    f = gvf.fabrique_cf(B0)
    f,g = gvf.fabrique_cf(B0, chi)
    q = gvf.fabrique_cf(B0, chi)



    import matplotlib.pyplot as plt
    plt.imshow(D.todense())
    plt.colorbar(orientation='horizontal')
    plt.title('D')
    plt.show()

    plt.plot(gvf.fabrique_q(B0))
    plt.title('q')
    plt.show()

    plt.imshow(M.todense())
    plt.colorbar()
    plt.title('M')
    plt.show()

    plt.imshow(G.todense())
    plt.colorbar()
    plt.title('G')
    plt.show()
