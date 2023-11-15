#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "methodesnum.h"
#include "donnees.h"

double k(double x, double t) {
    /* fonction qui calcule k(x,t) pour des reels x et t donnes */
    int M = recup_M(M);
    double H = recup_H(H), L = recup_L(L);
    int m;
    double K = 1./(H*L);
    for (m = 1; m<=M; m++) {
        K = K + (2*M_PI/pow(L,2))*m*cos(m*M_PI*x/L)*cos(m*M_PI*t/L)/sinh(m*M_PI*H/L);
    }
    return K;
}

double gauss_approx(double x, int n, double a, double b, double val[]) {
    /* fonction pour la quadrature de Gauss-Legendre dans les méthodes d'Adomain et des noyaux iteres
    x: point x pour lequel on calcule \int_a^b k(x,t) u_n(t) dt
    n: nombre de points de quadrature
    a: borne inferieure de l'integrale
    b: borne superieure de l'integrale
    val[]: tableaux de valeurs que l'on calcule dans la fonction noyaux_iter */

    int n = recup_n(n);
    double Inte;
    if (sizeof(val) != n) {
        printf("%s", "Mauvaise taille du tableau");
    }
    else {
        int i;
        float *T = malloc(sizeof(int[n])); // allocation dynamique du tableau T qui contient les points de quadrature
        float *Tp = malloc(sizeof(int[n])); // allocation dynamique du tableau T qui contient les poids de gauss
        switch (n) { // on remplit les tableaux pour les poids et points de quadrature en fonction de n
            case 1:
                T[0] = 0.;
                Tp[0] = 2.;
                break;
            case 2:
                T[0] = -1./sqrt(3.);
                T[1] = 1./sqrt(3.);
                Tp[0] = 1.;
                Tp[1] = 1.;
                break;
            case 3:
                T[0] = -sqrt(3.)/sqrt(5.);
                T[1] = 0.;
                T[2] = sqrt(3.)/sqrt(5.);
                Tp[0] = 5./9.;
                Tp[1] = 8./9.;
                Tp[2] = 5./9.;
                break;
            case 4:
                T[0] = -sqrt(3./7.-2./7.*sqrt(6./5.));
                T[1] = -sqrt(3./7.+2./7.*sqrt(6./5.));
                T[2] = sqrt(3./7.-2./7.*sqrt(6./5.));
                T[3] = sqrt(3./7.+2./7.*sqrt(6./5.));
                Tp[0] = (18.+sqrt(30.))/36;
                Tp[1] = (18.-sqrt(30.))/36;
                Tp[2] = (18.+sqrt(30.))/36;
                Tp[3] = (18.-sqrt(30.))/36;
                break;
            case 5:
                T[0] = 0;
                T[1] = -(1./3.)*sqrt(5.-2.*sqrt(10./7.));
                T[2] = -(1./3.)*sqrt(5.+2.*sqrt(10./7.));
                T[3] = (1./3.)*sqrt(5.-2.*sqrt(10./7.));
                T[4] = (1./3.)*sqrt(5.+2.*sqrt(10./7.));
                Tp[0] = 128./155.;
                Tp[1] = (322.+13.*sqrt(70.))/900.;
                Tp[2] = (322.-13.*sqrt(70.))/900.;
                Tp[3] = (322.+13.*sqrt(70.))/900.;
                Tp[4] = (322.-13.*sqrt(70.))/900.;
        }
        Inte = 0.;
        for (i = 0; i < n; i++) { // on calcule notre integrale
            Inte = Inte + ((b-a)/2.)*(Tp[i] * k(x,((b-a)*T[i]/2.)+(a+b)/2.)*val[i]);
        }
        // on desalloue la place memoire des tableaux
        free(T);
        free(Tp);
    }
    return Inte; // on retourne la valeur de l'integrale
}

/*double fctn_int(double x, double (*k)(double, double), double val[]) {
    int n = recup_n(n);
    if (sizeof(val) != n) {
        printf("%s", "Mauvaise taille du tableau");
    }
    else {

    }
}*/

double noyaux_iter(double x, double alpha) {
    double H = recup_H(H), L = recup_L(L);
    int n = recup_n(n); 
    

}