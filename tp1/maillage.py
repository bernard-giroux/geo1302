# -*- coding: utf-8 -*-
"""
Copyright 2016 Bernard Giroux

email: Bernard.Giroux@ete.inrs.ca

"""

import sys
import numpy as np

class TriSurf:
    """
    Classe pour représenter une surface en 3D avec un maillage triangulaire
    """
    def __init__(self, name=None):
        self.name = name
        self.p = []
        self.t = []

    def chargerMsh(self, fichier):
        """
        Méthode pour construire le maillage à partir d'un fichier produit avec
        le mailleur gmsh (http://gmsh.info)

        Le maillage doit avoir été sauvegardé en format 2 ASCII
        """
        try:
            with open(fichier, 'rt') as f:
                for line in f:
                    if line.find('$MeshFormat')>=0:
                        tmp = f.readline().split()
                        if tmp[0] != '2.2':
                            print('Erreur: format de fichier non reconnu')
                            sys.exit(1)
                        if tmp[1] != '0':
                            print('Erreur: fichier ASCII attendu')
                            sys.exit(1)

                        f.readline();  # $EndMeshFormat
                    elif line.find('$PhysicalNames')>=0:
                        nPhysicalNames = int( f.readline() )
                        for n in range(nPhysicalNames):
                            f.readline();
                        f.readline();  # $EndPhysicalNames
                    elif line.find('$Nodes')>=0:
                        nNodes = int( f.readline() )
                        self.p = np.empty([nNodes,3])
                        for n in range(nNodes):
                            tmp = f.readline().split()
                            self.p[n,0] = float( tmp[1] )
                            self.p[n,1] = float( tmp[2] )
                            self.p[n,2] = float( tmp[3] )
                        f.readline();  # $EndNodes
                    elif line.find('$Elements')>=0:
                        nElements = int( f.readline() )
                        self.t = np.empty([nElements,3],dtype=int)
                        ne = 0
                        for n in range(nElements):
                            tmp = f.readline().split()
                            if tmp[1] == '2':  # on ne charge que les triangles
                                nTags = int( tmp[2] )
                                self.t[ne,0] = int( tmp[3+nTags] )-1  # les indices commencent à 0 en python
                                self.t[ne,1] = int( tmp[4+nTags] )-1
                                self.t[ne,2] = int( tmp[5+nTags] )-1
                                ne += 1
                        self.t = self.t[:ne,:]  # on redimensionne en fct du nombre de triangles

            if len(self.p) == 0:
                print('Erreur: le fichier ne contient pas de noeuds')
                sys.exit(1)
            if len(self.t) == 0:
                print('Erreur: le fichier ne contient pas de triangles')
                sys.exit(1)

        except IOError:
            print("Erreur: impossible d'ouvrir le fichier "+fichier)
            sys.exit(1)


    def calculVolume(self,checkNormal=True):
        """
        Calcule le volume du polyèdre
        """
        V = 0.0
        Nf = self.t.shape[0]

        if checkNormal==True:
            x0 = np.sum(self.p[:,0])/self.p.shape[0]
            y0 = np.sum(self.p[:,1])/self.p.shape[0]
            z0 = np.sum(self.p[:,2])/self.p.shape[0]

        Un = np.empty((Nf,3))
        for t in range(Nf):
            v1 = self.p[self.t[t,2],:] - self.p[self.t[t,0],:]
            v2 = self.p[self.t[t,1],:] - self.p[self.t[t,0],:]
            ss = np.array([v2[1]*v1[2]-v2[2]*v1[1], v2[2]*v1[0]-v2[0]*v1[2], v2[0]*v1[1]-v2[1]*v1[0]])
            Un[t,:] = ss/np.linalg.norm(ss)

            if checkNormal==True:
                v1 = self.p[self.t[t,0],:] - np.array([x0, y0, z0])
                if np.sum(Un[t,:]*v1)<0:  # pointing toward centroid, change order
                    tmp=self.t[t,1]
                    self.t[t,1]=self.t[t,2]
                    self.t[t,2]=tmp
                    Un[t,:] = -Un[t,:]
            V += np.dot(self.p[self.t[t,0],:],np.cross(self.p[self.t[t,1],:],self.p[self.t[t,2],:]))

        return V/6.0

    def toVTK(self, fichier):
        """
        Sauvegarde en format VTK

        INPUT
            fichier: nom du fichier de sortie (string)
        """
        if len(self.p) == 0 or len(self.t) == 0:
            print("Erreur: maillage non défini")
            sys.exit(1)

        if fichier[-4:] != '.vtu':
            fichier = fichier + '.vtu'
            print(fichier)

        try:
            with open(fichier, 'wt') as f:
                f.write('<?xml version="1.0"?>\n')
                f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
                f.write('  <UnstructuredGrid>\n')
                f.write('    <Piece NumberOfPoints="{0:d}" NumberOfCells="{1:d}">\n'.format(self.p.shape[0], self.t.shape[0]))
                f.write('      <Points>\n')
                f.write('        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
                for n in np.arange(self.p.shape[0]):
                    f.write('          {0:f} {1:f} {2:f}\n'.format(self.p[n,0], self.p[n,1], self.p[n,2]))
                f.write('        </DataArray>\n')
                f.write('      </Points>\n')
                f.write('      <Cells>\n')
                f.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
                for n in np.arange(self.t.shape[0]):
                    f.write('          {0:d} {1:d} {2:d}\n'.format(self.t[n,0], self.t[n,1], self.t[n,2]))
                f.write('        </DataArray>\n')
                f.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
                off = 3
                for n in np.arange(self.t.shape[0]):
                    f.write('          {0:d}\n'.format(off))
                    off += 3
                f.write('        </DataArray>\n')
                f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
                for n in np.arange(self.t.shape[0]):
                    f.write('          5\n')
                f.write('        </DataArray>\n')
                f.write('      </Cells>\n')
                f.write('    </Piece>\n')
                f.write('  </UnstructuredGrid>\n')
                f.write('</VTKFile>\n')
        except IOError:
            print("Erreur: impossible d'ouvrir le fichier "+fichier)
            sys.exit(1)


if __name__ == '__main__':
    m = TriSurf()
    m.chargerMsh('/Users/giroux/Cours/modelisation/gmsh/sphere.msh')
    print(m.p.shape)
    print(m.t.shape)
    m.toVTK('sphere')
