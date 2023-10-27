#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "methodesnum.h"




double gauss(int n, double (*f)(double,int,double), double a, double b, int m, double alpha) {
    /* fonction pour calculer une integrale numeriquement par la methode de Gauss-Legendre
    n: nombre de points de quadrature: 1, 2, 3, 4 ou 5
    f: la fonction dont on souhaite approcher l'integrale
    a: borne inf du domaine sur lequel on integre
    b: borne sup du domaine sur lequel on integre
    retourne: la valeur approchee de l'integrale
    m: pour les sommes (on recupere l'indice m)
    alpha: parametre de Lavrentier qu'on recupere pour calculer f_3 */

    int i;
    double Inte;
    float *T = malloc(sizeof(int[n])); // allocation dynamique du tableau T qui contient les points de quadrature
    float *Tp = malloc(sizeof(int[n])); // allocation dynamique du tableau T qui contient les poids de gauss

    if (n!=1 && n!=2 && n!=3 && n!=4 && n!=5) { 
        printf("%s", "Nombre de points de quadrature non defini (n = 1, 2, 3, 4 ou 5)");
    }
    else {
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
            Inte = Inte + ((b-a)/2.)*(Tp[i] * f(((b-a)*T[i]/2.)+(a+b)/2.,m,alpha));
        }
    }

    // on desalloue la place memoire des tableaux
    free(T);
    free(Tp);

    return Inte; // on retourne la valeur de l'integrale
}

/*double newton(double x, double (*fonction)(double,double,double), double (*derivee)(double,double,double), double eps, double alpha) {
    /* fonction pour l'algorithme de Newton 
    x: point de Gamma_0 fixe pour lequel on cherche y
    fonction: fonction pour laquelle on cherche y telle que f(x,y) = 0
    derivee: derivee par rapport a y de notre fonction
    eps: tolerance epsilon pour l'approximation de y 

    double H = recup_H(H);
    double a = H/2.; // on commence au milieu
    double y;
    double delta = 1.;
    int N = 1000, i = 0;

    while (delta > eps && i<N) { // on donne un nombre max d'iterations
        y = -fonction(x,a,alpha)/derivee(x,a,alpha)+a;
        delta = fabs(y-a);
        a = y;
        i = i+1;
    }
    return a;
}*/