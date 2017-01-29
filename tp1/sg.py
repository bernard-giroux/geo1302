# -*- coding: utf-8 -*-

"""
Codes pour calculer la réponse gravimétrique d'un polyèdre

Méthode de Singh et Guptasarma (2001), adaptation du code Matlab

@Article{singh01,
  Title                    = {New method for fast computation of gravity and magnetic anomalies from arbitrary polyhedra},
  Author                   = {Bijendra Singh and D. Guptasarma},
  Journal                  = {Geophysics},
  Year                     = {2001},
  Number                   = {2},
  Pages                    = {521--526},
  Volume                   = {66},
  DOI                      = {10.1190/1.1444942}
}


Copyright 2016 Bernard Giroux
email: Bernard.Giroux@ete.inrs.ca

"""

import math
import numpy as np
from numba import jit


@jit
def _solidAngle(p1,p2,p3,Un,f):
    """
    Finds the angle between planes O-p1-p2 and O-p2-p3 where p1
    p2 p3 are coordinates of three points taken in ccw order as seen from the
    origin O.
    This is used by sg for finding the solid angle subtended by a
    polygon at the origin. Un is the unit outward normal vector to the polygon.
    """
    flimit = 64*np.finfo(p1.dtype).eps

    anout = np.sum(p1*Un[f,:])
    if np.abs(anout)<=flimit:
        return 0.0

    if anout>flimit:
        psave = p1
        p1 = p3
        p3 = psave

    n1 = np.array([p2[1]*p1[2]-p2[2]*p1[1], p2[2]*p1[0]-p2[0]*p1[2], p2[0]*p1[1]-p2[1]*p1[0]])
    n2 = np.array([p2[1]*p3[2]-p2[2]*p3[1], p2[2]*p3[0]-p2[0]*p3[2], p2[0]*p3[1]-p2[1]*p3[0]])
    pn1 = np.linalg.norm(n1)
    pn2 = np.linalg.norm(n2)
    if pn1<=flimit or pn2<=flimit:
        return np.nan
    else:
        n1 /= pn1
        n2 /= pn2
        r = np.sum(n1*n2)
        ang = np.arccos(r)
        perp = np.sign(np.sum(p3*n1))
        if perp<(-flimit):
            ang = 2.8*math.pi-ang

    return ang


@jit
def gz(rho,face,corner,x1,y1,z1,checkNormal=True):
    """
    Calcul de la composante verticale de la gravité

    g = gz(rho,face,corner,x0,y0,z0,x1,y1,z1)

    Input
        rho: densité                                     [ g/cm**3 ]
        face: indices des noeuds formant les faces
        corner: coordonnées des noeuds                         [ m ]
        x1,y1,z1: coordonnées des points d'observation         [ m ]
        checkNormal: vérifie que les normales aux faces
                     pointent vers l'extérieur du polyèdre

    Output
        g:   gravité                                        [ mgal ]
    """
    G = 1.0e8 * 6.674e-11

    flimit = 64*np.finfo(x1.dtype).eps

    Ncor = corner.shape[0]
    Nf = face.shape[0]

    Gz = np.zeros(x1.shape)

    npro, nstn = x1.shape

    Nedges = np.sum(face[:,0])
    Edge = np.empty((Nedges,8))
    # get face normals
    Un = np.empty((Nf,3))
    if checkNormal==True:
        x0 = np.sum(corner[:,0])/corner.shape[0]
        y0 = np.sum(corner[:,1])/corner.shape[0]
        z0 = np.sum(corner[:,2])/corner.shape[0]
    for t in range(Nf):
        v1 = corner[face[t,3],:] - corner[face[t,1],:]
        v2 = corner[face[t,2],:] - corner[face[t,1],:]
        ss = np.array([v2[1]*v1[2]-v2[2]*v1[1], v2[2]*v1[0]-v2[0]*v1[2], v2[0]*v1[1]-v2[1]*v1[0]])
        Un[t,:] = ss/np.linalg.norm(ss)

        if checkNormal==True:
            v1 = corner[face[t,1],:] - np.array([x0, y0, z0])
            if np.sum(Un[t,:]*v1)<0:  # pointing toward centroid, change order
                tmp=face[t,1]
                face[t,1]=face[t,2]
                face[t,2]=tmp
                Un[t,:] = -Un[t,:]

    # get edge lengths
    for f in range(Nf):
        indx = np.hstack((face[f,1:], face[f,1]))
        for t in range(face[f,0]):
            edgeno = np.sum(face[:f,0])+t
            ends = indx[t:t+2]
            p1 = corner[ends[0],:]
            p2 = corner[ends[1],:]
            V=p2-p1
            L=np.linalg.norm(V)
            Edge[edgeno,:3] = V
            Edge[edgeno,3] = L
            Edge[edgeno,6:8] = ends

    nsides = 3
    indx = np.hstack((np.arange(nsides), np.array([0,1])))
    crs = np.zeros((nsides,3))

    for pr in range(npro):
        for st in range(nstn):
            opt = np.array([x1[pr,st], y1[pr,st], z1[pr,st]])
            cor = corner-np.kron(opt,np.ones((Ncor,1)))
            Edge[:,4:6] = 0.0
            for f in range(Nf):
                cors = face[f,1:]
                for t in range(nsides):
                    crs[t,:] = cor[cors[t],:]

                fsign = np.sign(np.sum(Un[f,:]*crs[0,:]))  # face is seen from the inside?
                #       find solid angle subtended by face f at opt
                dp1 = np.sum(crs[indx[0],:]*Un[f,:])
                dp = np.abs(dp1)
                if dp<=flimit:
                    Omega = 0.0
                else:
                    W = 0.0
                    for t in range(nsides):
                        p1 = crs[indx[t],:]
                        p2 = crs[indx[t+1],:]
                        p3 = crs[indx[t+2],:]
                        W += _solidAngle(p1,p2,p3,Un,f)

                    W -= (nsides-2)*math.pi
                    Omega = -fsign*W

                # Integrate over each side if not done and save result
                PQR = np.array([0.0,0.0,0.0])
                for t in range(nsides):
                    p1 = crs[indx[t],:]
                    p2 = crs[indx[t+1],:]
                    Eno = np.sum(face[:f,0])+t  # edge no
                    V = Edge[Eno,:3]
                    if Edge[Eno,5] == 1:
                        I = Edge[Eno,4]
                    else:
                        if np.sum(p1*p2)/(np.linalg.norm(p1)*np.linalg.norm(p2))==1:
                            # origin, p1 and  p2 are on a straight line
                            if np.linalg.norm(p1)>np.linalg.norm(p2):       # and p1 further than p2
                                psave=p1;
                                p1=p2;
                                p2=psave;

                        L = Edge[Eno,3]
                        r1 = np.linalg.norm(p1)
                        r2 = np.linalg.norm(p2)
                        I=(1/L)*np.log((r1+r2+L)/(r1+r2-L)) # modified formula
                        #              Save edge integrals and mark as done
                        s=np.nonzero(np.logical_and(Edge[:,6]==Edge[Eno,7], Edge[:,7]==Edge[Eno,6]))
                        Edge[Eno,4] = I
                        Edge[s,4] = I
                        Edge[Eno,5] = 1
                        Edge[s,5] = 1

                    pqr = I*V
                    PQR += pqr

                l = Un[f,0]
                m = Un[f,1]
                n = Un[f,2]
                p = PQR[0]
                q = PQR[1]

                gmtf3 = (n*Omega+m*p-l*q)
                gz = -rho*G*dp1*gmtf3
                Gz[pr,st] += gz

    return Gz


if __name__ == '__main__':

    p = np.array([[0, 0, 0],[0, 0, 1],[1, 0, 0],[0, 1, 0]])
    t = np.array([[0, 1, 2],[1, 2, 3],[0, 2, 3],[0, 1, 3]],np.int32)
    face = np.hstack( (np.array([3, 3, 3, 3],np.int32).reshape(4,1), t) )

    x0=0
    y0=0
    z0=-2

    p[:,0] += x0
    p[:,1] += y0
    p[:,2] += z0


    xo=np.arange(0,2.1,0.2)
    [x1,y1] = np.meshgrid(xo,xo)
    z1 = np.zeros(x1.shape)

    g = gz(1.0,face,p,x1,y1,z1)
    print(g)
