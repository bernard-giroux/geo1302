# -*- coding: utf-8 -*-
import time
import numpy as np
import h5py

def prd(rho, x0, x, y, z):
    """
    PRD - Réponse gravimétrique d'un prisme rectangulaire droit

    g = prd(rho, x0, x, y, z)

    Input
        rho: densité                                     [ g/cm**3 ]
        x0:  coordonnées [x y z] du point d'observation        [ m ]
        x:   coord inférieure et supérieure du prisme selon x  [ m ]
        y:   coord inférieure et supérieure du prisme selon y  [ m ]
        z:   coord inférieure et supérieure du prisme selon z  [ m ]

    Output
        g:   gravité                                        [ mgal ]
    """
    G = 1.0e8 * 6.674e-11
    g = 0.0
    for i in [0, 1]:
        xi = x0[0]-x[i]
        for j in [0, 1]:
            yj = x0[1]-y[j]
            for k in [0, 1]:
                zk = x0[2]-z[k]
                mu =
                rijk =
                g +=

    return -G*rho*g

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
        return (

    def prd_G(self, x0):
        """
        PRD_G - Opérateur direct gravimétrique pour une grille de prismes rectangulaires droits

        G = prd_G(x0)

        Input
            x0: coordonnées des points d'observation (N x 3)

        Output
            G: opérateur direct (N x M)
        """
        M = self.nc
        N = x0.shape[0]
        G = np.ndarray((N,M))

        for k in np.arange(self.nz):
            for j in np.arange(self.ny):
                for i in np.arange(self.nx):
                    m = self.ind(i, j, k)
                    for n in np.arange(N):
                        G[n,m] = prd(
        return G

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




if __name__ == '__main__':

    # test de la fonction prd
    rho = 0.2
    x = (10, 15)
    y = (20, 25)
    z = (5, 15)

    x0 = (0, 0, 0)

    g = prd(rho, x0, x, y, z)
    print(g)

    x0 = (12.5, 22.5, 10)

    g = prd(rho, x0, x, y, z)
    print(g)

    x = np.arange(-8.5, 9.0)
    y = np.arange(-10.5, 11.0)
    z = np.arange(10.0)
    g = Grille(x, y, z)

    print(g.nc, g.dx)

    print(g.ind(0, 0, 0))
    print(g.ind(1, 0, 0))

    rho = np.zeros((g.nc,))
    rho[g.ind(1,1,1)] = 1.0
    rho[g.ind(1,2,1)] = 2.0
    rho[g.ind(1,1,3)] = 3.0

    x0 = np.array([[0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [5.0, 0.0, 0.0]])

    tic = time.time()
    G = g.prd_G(x0)
    t_G = time.time() - tic

    rho = np.zeros((g.nc,))
    rho[ g.ind(8,10,5) ] = 10.0

    tic = time.time()
    gz = np.dot(G, rho)
    t_mult = time.time() - tic

    print(t_G, t_mult)

    import matplotlib.pyplot as plt
    plt.plot(x0[:,0], gz, 'o')
    plt.xlabel('x (m)', fontsize=12)
    plt.ylabel('$g_z$ (mgal)', fontsize=12)

    g.toXdmf(rho, 'rho', 'grille')
