#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fonctions.h"
#include "methodesnum.h"
#include "donnees.h"

double alpha = 1.;
double L = recup_L(L);

/* Discrétisation de [0,L]*/

int D = 10; /* Nombre de points sur la discrétisation de [0,L] */
double dx = L/D;

double *xl = malloc(sizeof(int[D]));



double u0(double x){
	double res;
	res = (1./alpha)*h(x);
	return res;
}

/* Création du tableau U1 qui contient l'évaluation de u_{n+1} en x_l */

double *U1 = malloc(sizeof(int[D]));

/* Création du tableau U qui contient l'évaluation de u_n en x_l */

double *U = malloc(sizeof(int[D]));

/* Initialisation de U  */

for(l=0;l<D+1;l++){
	U[l]=u0(xl[l]);
}

/* Le fichier doit renvoyer la fonction f3 app par la méthode d'Adomain sous le format double f3_Ad(x) */

/* Créer une fonction pour intégrer la fonction calculée par la méthodes de Gauss */

