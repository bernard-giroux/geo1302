
# coding: utf-8

# # Modélisation de la température en régime permanent

import matplotlib.pyplot as plt
import numpy as np


# ## Définition des paramètres

T0 = 280      # température en surface (K)
Q = -55e-3    # flux dans le N^e couche (W/m^2)
N = 20        # Nombre de couches
dz = 100      # épaisseur des couches (m)

# Production de chaleur
A = np.zeros((N,1))
A[-3:-1] = 2.8e-6       # W/m^3

# Conductivité thermique  (W/m/K)
l = np.zeros((N,1))
l[:5] = 1.8
l[5:10] = 3.7
l[10:16] = 2.4
l[-4:] = 3.5


y=np.vstack([np.array([0]),np.kron((np.arange(1,N)[:,np.newaxis]*dz),
                                   np.array([[1],[1]])),np.array([N*dz])])

fig = plt.figure()
plt.subplot(121)
plt.plot(np.kron(l,np.array([[1],[1]])),y);
plt.gca().invert_yaxis()
plt.xlabel(r'$\lambda$ (W/m/K)',fontsize=14);
plt.ylabel('Profondeur (m)',fontsize=14);

plt.subplot(122)
plt.plot(np.kron(1000*A,np.array([[1],[1]])),y);
plt.gca().invert_yaxis()
plt.xlabel('A (mW/m$^3$)',fontsize=14);
#plt.ylabel('Profondeur (m)',fontsize=14);
plt.xlim([-0.001, 0.003])
plt.xticks(np.array([0, 2e-3]));
plt.show()

# ## Construction de la matrice A

M = np.zeros(())
M[0,0] = 1
for n in np.arange(1,N):
    M[n,n-1] =
    M[n,n] =
    M[n,n+1] =

M[N,N-1] = l[N-1]/dz
M[N,N] = -l[N-1]/dz

plt.spy(M,marker='o');
plt.show()

# ### Affichage de la matrice

m=np.unique(M)
MM = M.copy()
nn = np.arange(m.size)
label = []
for n in nn:
    MM[M==m[n]] = n
    label.append(str(m[n]))

plt.matshow(MM);
h = plt.colorbar(ticks=nn);
h.ax.set_yticklabels(label);
plt.show()

# ## Construction du vecteur source

b = np.zeros((N+1,1))
b[0] =
for n in np.arange(1,N):
    b[n] = -(A[n-1]+A[n])/2

b[-1] =


# ## Solution du système

x = np.linalg.solve(M,b)

z = np.arange(0,N+1)[:,np.newaxis]*dz

plt.plot(x,z)
plt.gca().invert_yaxis()
plt.xlabel('T (K)',fontsize=14);
plt.ylabel('Profondeur (m)',fontsize=14);
plt.show()

y = b-M.dot(x)

plt.plot(z,y,'o-');
plt.xlabel('Profondeur (m)',fontsize=14);
plt.ylabel('Erreur',fontsize=14);
plt.show()

help(np.linalg.solve)
