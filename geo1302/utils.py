# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 19:54:19 2017

@author: giroux
"""

import numpy as np
import h5py
import inspect,dis

def nargout():
    """
    Return how many values the caller is expecting

    taken from
    http://stackoverflow.com/questions/16488872/python-check-for-the-number-of-output-arguments-a-function-is-called-with
    """
    f = inspect.currentframe()
    f = f.f_back.f_back
    c = f.f_code
    i = f.f_lasti
    bytecode = c.co_code
    instruction = bytecode[i+3]
    if instruction == dis.opmap['UNPACK_SEQUENCE']:
        howmany = bytecode[i+4]
        return howmany
    elif instruction == dis.opmap['POP_TOP']:
        return 0
    return 1

class Grille:
    """
    Classe pour gérer des grilles régulières (prismes rectangulaires droits)
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
        Retourne l'indice d'un prisme

        Input
            i: indice du voxel selon x
            j: indice du voxel selon y
            k: indice du voxel selon z
        """
        return (k * self.ny + j) * self.nx + i

    def __update_n(self):
        if '_x' in self.__dict__ and '_y' in self.__dict__ and '_z' in self.__dict__:
            self.nx = len(self._x)-1             # nombre de prismes selon x
            self.ny = len(self._y)-1             # nombre de prismes selon y
            self.nz = len(self._z)-1             # nombre de prismes selon z
            self.nc = self.nx * self.ny * self.nz    # nombre de prismes

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
            if len(np.unique(np.diff(tmp)))>1:
                raise ValueError('Constant step size needed')
        except:
            raise

        self._x = tmp
        self.dx = tmp[1]-tmp[0]
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
            if len(np.unique(np.diff(tmp)))>1:
                raise ValueError('Constant step size needed')
        except:
            raise

        self._y = tmp
        self.dy = tmp[1]-tmp[0]
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
            if len(np.unique(np.diff(tmp)))>1:
                raise ValueError('Constant step size needed')
        except:
            raise

        self._z = tmp
        self.dz = tmp[1]-tmp[0]
        self.__update_n()

    def toXdmf(self, field, fieldname, filename):
        """
        Save a field in xdmf format (http://www.xdmf.org/index.php/Main_Page)

        INPUT
            field: data array of size equal to the number of cells in the grid
            fieldname: name to be assinged to the data (string)
            filename: name of xdmf file (string)
        """
        ox = self.x[0]# + self.dx/2
        oy = self.y[0]# + self.dy/2
        oz = self.z[0]# + self.dz/2

        f = open(filename+'.xmf','w')

        f.write('<?xml version="1.0" ?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        f.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
        f.write(' <Domain>\n')
        f.write('   <Grid Name="Structured Grid" GridType="Uniform">\n')
        f.write('     <Topology TopologyType="3DCORECTMesh" NumberOfElements="'+
        repr(self.nz+1)+' '+repr(self.ny+1)+' '+repr(self.nx+1)+'"/>\n')
        f.write('     <Geometry GeometryType="ORIGIN_DXDYDZ">\n')
        f.write('       <DataItem Dimensions="3 " NumberType="Float" Precision="4" Format="XML">\n')
        f.write('          '+repr(oz)+' '+repr(oy)+' '+repr(ox)+'\n')
        f.write('       </DataItem>\n')
        f.write('       <DataItem Dimensions="3 " NumberType="Float" Precision="4" Format="XML">\n')
        f.write('        '+repr(self.dz)+' '+repr(self.dy)+' '+repr(self.dx)+'\n')
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




if __name__ == '__main__':

    x = [1, 2, 3, 4]
    y = [0, 1, 2, 3, 4, 5, 6, 7]
    z = np.arange(-5,5)
    g = Grille(x, y, z)

    print(g.nc, g.dx)

    print(g.ind(1,1,1))

    rho = np.zeros((g.nc,))
    rho[g.ind(1,1,1)] = 1.0
    rho[g.ind(1,2,1)] = 2.0
    rho[g.ind(1,1,3)] = 3.0

    g.toXdmf(rho, 'rho', 'rho')
