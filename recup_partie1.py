# recuperation des donnees pour le projet 2 - calcul numerique scientifique
# Partie 1 -resolution du probleme de Cauchy

#%% importation packages

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import axes3d
import os

#%% definition des variables

L = 1.
H = 1.

Alpha = np.zeros(22)

Alpha[0] = 10.
Alpha[1] = 1.
for i in range(2,22):
    Alpha[i] = 10**(-(i-1))
    
N = 51 #nombre de points de maillage pour x
    
X = np.linspace(0.,L,N)
Y = np.linspace(0.,H,N)
  
#%% execution du programme partie1

#os.system('make clean')
#os.system('make partie1')
#os.system('./partie1') # on execute le fichier executable partie1 depuis notre console python
# penser a le compiler d'abord
    

#%% lecture des fichiers

# solution exacte sur Gamma_3
file = open('exact_1.txt', 'r')

data = file.read()
Exact = data.split()
Ex = [float(i) for i in Exact]

#print(Ex)


# solution approchee sur Gamma_3
file = open('approche_noyaux.txt', 'r')

data = file.read()
Approche = data.split()
Approche = [float(i) for i in Approche]
App = np.array(Approche)

#print(App)

# On divise le tableau en 20 tableaux (1 tableau pour chaque alpha)
T_alpha = np.split(App, 22)

# solution approchee sur Gamma_3
file = open('approche_adomain.txt', 'r')

data = file.read()
Approche = data.split()
Approche = [float(i) for i in Approche]
App_ado = np.array(Approche)

#print(App)

# On divise le tableau en 20 tableaux (1 tableau pour chaque alpha)
T_alpha_ado = np.split(App_ado, 22)

#%% calcul des erreurs

Err = np.zeros(22)
for i in range(22):
    Err[i] = np.linalg.norm(T_alpha[i] - Ex)

Err_ado = np.zeros(22)
for i in range(22):
    Err_ado[i] = np.linalg.norm(T_alpha_ado[i] - Ex)/ np.linalg.norm(Ex)

#%% graphiques

plt.loglog(Alpha, Err)
plt.title('Erreur entre les solutions exacte et approchée (noyaux) sur Gamma_3')
plt.show()

plt.loglog(Alpha, Err_ado)
plt.title('Erreur entre les solutions exacte et approchée (Adomain) sur Gamma_3')
plt.show()

plt.plot(X, Ex)
plt.title('Solution exacte T_ex sur Gamma_3')
plt.show()

plt.plot(X, T_alpha[1])
plt.title('Solution approchee par la méthode des noyaux itérés T_app sur Gamma_3 pour alpha = 1.0')
plt.show()

plt.plot(X, T_alpha_ado[1])
plt.title('Solution approchee par la méthode d Adomain T_app sur Gamma_3 pour alpha = 1.0')
plt.show()


#%% 

def f(x,y):
    return np.cosh(np.pi*y)*np.cos(np.pi*x)

#plt.plot(X, f_3(X))
#plt.show()


#%% lecture des fichiers pour un alpha donne = 1.

alpha = 1.

file = open('solapp_omega_noyaux.txt', 'r')
data = file.read()
SolApp = data.split()
SolApp = [float(i) for i in SolApp]
SolApp = np.array(SolApp)
SolApp = np.split(SolApp, N)
Z = np.vstack(SolApp)

file = open('solapp_omega_adomain.txt', 'r')
data = file.read()
SolApp = data.split()
SolApp = [float(i) for i in SolApp]
SolApp = np.array(SolApp)
SolApp = np.split(SolApp, N)
Z_ado = np.vstack(SolApp)

file = open('solex_omega.txt', 'r')
data = file.read()
SolEx = data.split()
SolEx = [float(i) for i in SolEx]
SolEx = np.array(SolEx)
SolEx = np.split(SolEx, N)
Zex = np.vstack(SolEx)

file = open('erreur_omega_noyaux.txt', 'r')
data = file.read()
ERR = data.split()
ERR = [float(i) for i in ERR]
ERR = np.array(ERR)
ERR = np.split(ERR, N)
Zerr = np.vstack(ERR)

file = open('erreur_omega_adomain.txt', 'r')
data = file.read()
ERR = data.split()
ERR = [float(i) for i in ERR]
ERR = np.array(ERR)
ERR = np.split(ERR, N)
Zerr_ado = np.vstack(ERR)


#%% graphiques pour un alpha donne

# Tracé du résultat en 3D
#ax = plt.figure().add_subplot(projection = '3d')
#X,Y = np.meshgrid(X,Y)
#ax.plot_surface(X,Y,Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')


# solution approchee
clev = np.arange(Z.min(),Z.max(),0.1)
h = plt.contourf(X, Y, Z, clev)
plt.axis('scaled')
plt.colorbar()
plt.title('Solution approchee (noyaux) pour alpha = '+str(alpha))
plt.xlabel('x')
plt.ylabel('y')
plt.show()

clev = np.arange(Z_ado.min(),Z_ado.max(),0.1)
h = plt.contourf(X, Y, Z_ado, clev)
plt.axis('scaled')
plt.colorbar()
plt.title('Solution approchee (Adomain) pour alpha = '+str(alpha))
plt.xlabel('x')
plt.ylabel('y')
plt.show()

#solution exacte
clev = np.arange(Zex.min(),Zex.max(),0.1)
h = plt.contourf(X, Y, Zex, clev)
plt.axis('scaled')
plt.colorbar()
plt.title('Solution exacte')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

#erreur
clev = np.arange(Zerr.min(),Zerr.max(),0.1)
h = plt.contourf(X, Y, Zerr, clev)
plt.axis('scaled')
plt.colorbar()
plt.title('Erreur: T_ex-T_cal (noyaux)')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

clev = np.arange(Zerr_ado.min(),Zerr.max(),0.1)
h = plt.contourf(X, Y, Zerr_ado, clev)
plt.axis('scaled')
plt.colorbar()
plt.title('Erreur: T_ex-T_cal (Adomain)')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

#X,Y = np.meshgrid(X,Y)
#h = plt.contourf(X, Y, f(X,Y))
#lt.axis('scaled')
#plt.colorbar()
#plt.title('Test python')
#plt.show()