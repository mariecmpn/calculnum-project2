#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "fonctions.h"
#include "methodesnum.h"
#include "donnees.h"
#include "noyaux.h"


double Adomain(double x, double alpha){

	/* On récupère les données du problème */
	
	double L = recup_L(L);
	int n = recup_n(n);
	int nb_iter = recup_nb_iter(nb_iter);
	
	/* Notre résultat est stocké dans res*/
	double res = 0.;
	
	/* On ne calcule qu'une seule fois 1/alpha */
	double a = 1./alpha;

    	int i,j; //entiers pour les boucles for

    	double u_n; // resultat de l'equation de fredholm de seconde espece

    	// on definit u_0
    	double *T = malloc(sizeof(int[n])); // allocation dynamique du tableau T qui contient les points de quadrature
    	double *U = malloc(sizeof(int[n])); // allocation dynamique du tableau U qui contiendra les valeurs de u_n-1 a chaque iteration n
    	double *Un = malloc(sizeof(int[n])); // allocation dynamique du tableau U qui contiendra les valeurs de u_n a chaque iteration n
    	
    	switch (n) { // on remplit les tableaux des points de quadrature, et de u_0 en fonction de n
    	case 1:
            	T[0] = 0.;
            	U[0] = a*h(T[0]);
            	break;
        case 2:
            	T[0] = -1./sqrt(3.);
            	T[1] = 1./sqrt(3.);
            	U[0] = a*h(T[0]);
            	U[1] = a*h(T[1]);
            	break;
        case 3:
            	T[0] = -sqrt(3.)/sqrt(5.);
            	T[1] = 0.;
            	T[2] = sqrt(3.)/sqrt(5.);
            	U[0] = a*h(T[0]);
            	U[1] = a*h(T[1]);
            	U[2] = a*h(T[2]);
            	break;
        case 4:
            	T[0] = -sqrt(3./7.-2./7.*sqrt(6./5.));
            	T[1] = -sqrt(3./7.+2./7.*sqrt(6./5.));
            	T[2] = sqrt(3./7.-2./7.*sqrt(6./5.));
            	T[3] = sqrt(3./7.+2./7.*sqrt(6./5.));
            	U[0] = a*h(T[0]);
            	U[1] = a*h(T[1]);
            	U[2] = a*h(T[2]);
            	U[3] = a*h(T[3]);
            	break;
        case 5:
            	T[0] = 0;
            	T[1] = -(1./3.)*sqrt(5.-2.*sqrt(10./7.));
            	T[2] = -(1./3.)*sqrt(5.+2.*sqrt(10./7.));
            	T[3] = (1./3.)*sqrt(5.-2.*sqrt(10./7.));
            	T[4] = (1./3.)*sqrt(5.+2.*sqrt(10./7.));
            	U[0] = a*h(T[0]);
            	U[1] = a*h(T[1]);
            	U[2] = a*h(T[2]);
            	U[3] = a*h(T[3]);
            	U[4] = a*h(T[4]);
        }

	for (i = 1; i <=nb_iter; i++) {
                u_n = - a*gauss_approx(x,n,0,L,U); // on calcule U_n au point x
                res = res + u_n;
            	for (j = 0; j<n; j++) { // puis on calcule u_n aux points de quadrature pour l'iteration d'apres
                	Un[j] = - a*gauss_approx(T[j],n,0,L,U);
            	}
                for (j = 0; j<n; j++) { 
                    U[j] = Un[j];
                }
    	}

    // on libere l'espace memoire des tableaux alloues dynamiquement
    
	free(T);
	free(U);
	free(Un);
	
	return res;
}

/* Le fichier doit renvoyer la fonction f3 app par la méthode d'Adomain sous le format double f3_Ad(x) */

/* Créer une fonction pour intégrer la fonction calculée par la méthodes de Gauss */


/*double adomain_direct(double x, double alpha) {
    int nb_iter = recup_nb_iter(nb_iter);
    int j; // entier pour boucle for
    double L = recup_L(L);
	int n = recup_n(n);
    /* On ne calcule qu'une seule fois 1/alpha 
	double a = 1./alpha;

    // on definit u_0
    	double *T = malloc(sizeof(int[n])); // allocation dynamique du tableau T qui contient les points de quadrature
    	double *U = malloc(sizeof(int[n])); // allocation dynamique du tableau U qui contiendra les valeurs de u_n-1 a chaque iteration n
    	double *Un = malloc(sizeof(int[n])); // allocation dynamique du tableau U qui contiendra les valeurs de u_n a chaque iteration n
    	
    	switch (n) { // on remplit les tableaux des points de quadrature, et de u_0 en fonction de n
    	case 1:
            	T[0] = 0.;
            	U[0] = a*h(T[0]);
            	break;
        case 2:
            	T[0] = -1./sqrt(3.);
            	T[1] = 1./sqrt(3.);
            	U[0] = a*h(T[0]);
            	U[1] = a*h(T[1]);
            	break;
        case 3:
            	T[0] = -sqrt(3.)/sqrt(5.);
            	T[1] = 0.;
            	T[2] = sqrt(3.)/sqrt(5.);
            	U[0] = a*h(T[0]);
            	U[1] = a*h(T[1]);
            	U[2] = a*h(T[2]);
            	break;
        case 4:
            	T[0] = -sqrt(3./7.-2./7.*sqrt(6./5.));
            	T[1] = -sqrt(3./7.+2./7.*sqrt(6./5.));
            	T[2] = sqrt(3./7.-2./7.*sqrt(6./5.));
            	T[3] = sqrt(3./7.+2./7.*sqrt(6./5.));
            	U[0] = a*h(T[0]);
            	U[1] = a*h(T[1]);
            	U[2] = a*h(T[2]);
            	U[3] = a*h(T[3]);
            	break;
        case 5:
            	T[0] = 0;
            	T[1] = -(1./3.)*sqrt(5.-2.*sqrt(10./7.));
            	T[2] = -(1./3.)*sqrt(5.+2.*sqrt(10./7.));
            	T[3] = (1./3.)*sqrt(5.-2.*sqrt(10./7.));
            	T[4] = (1./3.)*sqrt(5.+2.*sqrt(10./7.));
            	U[0] = a*h(T[0]);
            	U[1] = a*h(T[1]);
            	U[2] = a*h(T[2]);
            	U[3] = a*h(T[3]);
            	U[4] = a*h(T[4]);
        }

        for (i = 1; i <=nb_iter; i++) {
            	for (j = 0; j<n; j++) {
                	Un[j] = - a*gauss_approx(T[j],n,0,L,U);
            	}
                for (j = 0; j<n; j++) {
                    U[j] = Un[j];
                }
    	}

    // on libere l'espace memoire des tableaux alloues dynamiquement
    
	free(T);
	free(U);
	free(Un);
	
	return res;
}*/