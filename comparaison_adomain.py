# recuperation des donnees pour le projet 1 - calcul numerique scientifique
# Comparaison methode directe et recurrence pour Adomain

#%% importation packages

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import axes3d
import os

#%% definition des variables

L = 1.
H = 1.

    
N = 51 #nombre de points de maillage pour x
    
X = np.linspace(0.,L,N)

#%% recup donnees

file = open('adomain_direct.txt', 'r')

data = file.read()
Direct = data.split()
Dir = [float(i) for i in Direct]

file = open('adomain_indirect.txt', 'r')

data = file.read()
Indirect = data.split()
Indir = [float(i) for i in Indirect]

file = open('adomain_erreur.txt', 'r')

data = file.read()
Erreur = data.split()
Err = [float(i) for i in Erreur]

#%% graphiques

plt.plot(X, Dir, label = 'direct')
plt.plot(X, Indir, label = 'indirect')
plt.legend()
plt.title('Comparaison methode directe et recurrence d Adomain')
plt.show()

plt.plot(X, Err)
plt.title('Erreurs entre methode directe et recurrence pour chaque point de maillage')
plt.show()