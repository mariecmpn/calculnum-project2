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
    double l = 1./L; // au lieu de diviser a chaque fois on calcule 1/L une seule fois
    double K = 1./(H*L);
    for (m = 1; m<=M; m++) {
        K = K + (2*M_PI*l*l)*m*cos(m*M_PI*x*l)*cos(m*M_PI*t*l)/sinh(m*M_PI*H*l);
    }
    return K;
}


double h(double x) {
    /* fonction qui approche (par une somme de M termes) la fonction h
    x: reel dont on veut calculer l'image */
    double S;
    int i;
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    double l = 1./L; // au lieu de diviser a chaque fois on calcule 1/L une seule fois
    // puis on retourne la fonction souhaitee
    /*S = -q_0(x) + (1/(L*H))*gauss(n,f_0,0,L,0,0);
    for (i = 1; i<=M; i++) {
        //S = S + (2*M_PI/pow(L,2)) * (i*cos(i*M_PI*x/L)/sinh(i*M_PI*x/L)) * f_0(x,i,0.) * cosh(i*M_PI*H/L) * gauss(n,inte_h,0,L,i,0);
        S = S + (2*M_PI/pow(L,2)) * ((i*cos(i*M_PI*x/L)*cosh(i*M_PI*H/L))/sinh(i*M_PI*H/L))* gauss(n,inte_h,0.,L,i,0.);
    }*/
    S = (M_PI*l)*(cos(M_PI*x*l)*cosh(M_PI*H*l))/sinh(M_PI*H*l); // on enleve la somme pour cet exemple
    return S;
}

double inte_h(double x, int m, double alpha) {
    /* fonction a integrer dans l'expression de h */
    //on recupere d'abord les donnees du probleme
    double L,H;
    L = recup_L(L);
    H = recup_H(H);
    // puis on retourne la fonction souhaitee
    //double r = f_0(x, m, 0.)*cos(m*M_PI*x/L);
    double r = h(x);
    return r;
}

double gauss_approx(double x, int n, double a, double b, double val[5]) {
    /* fonction pour la quadrature de Gauss-Legendre dans les méthodes d'Adomain et des noyaux iteres
    x: point x pour lequel on calcule \int_a^b k(x,t) u_n(t) dt
    n: nombre de points de quadrature
    a: borne inferieure de l'integrale
    b: borne superieure de l'integrale
    val[]: tableaux de valeurs que l'on calcule dans la fonction noyaux_iter */

    double Inte;
    int i;
    float T[5]; // allocation dynamique du tableau T qui contient les points de quadrature
    float Tp[5]; // allocation dynamique du tableau T qui contient les poids de gauss
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
        Inte = Inte + ((b-a)*0.5)*(Tp[i] * k(x,((b-a)*T[i]*0.5)+(a+b)*0.5)*val[i]);
    }

    //free(val);
    return Inte; // on retourne la valeur de l'integrale
}

double noyaux_iter(double x, double alpha) {
    double H = recup_H(H), L = recup_L(L);
    int n = recup_n(n); 
    int i,j; //entiers pour les boucles for
    int nb_iter = recup_nb_iter(nb_iter);
    double a = 1./alpha; // on calcule le coefficient 1/alpha une seule fois
    double u_n; // resultat de l'equation de fredholm de seconde espece

    // on definit u_0
    double *T = malloc(sizeof(int[n])); // allocation dynamique du tableau T qui contient les points de quadrature
    double *U = malloc(sizeof(int[n])); // allocation dynamique du tableau U qui contiendra les valeurs de u_n-1 a chaque iteration n
    double *Un = malloc(sizeof(int[n])); // allocation dynamique du tableau U qui contiendra les valeurs de u_n a chaque iteration n
    switch (n) { // on remplit les tableaux des points de quadrature, et de u_0 en fonction de n
    // pour la premiere iteration on prend u_0(x) = 1/alpha * h(x)
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
        if (i == nb_iter) { // pour notre derniere iteration on calcule u_{n+1}
            u_n = a*h(x) - a*gauss_approx(x,n,0,L,U);
        } 
        else { // si ce n'est pas la derniere iteration on calcule les u_n(t) pour t les points de quadrature de gauss-legendre
            for (j = 0; j<n; j++) {
                Un[j] = a*h(T[j]) - a*gauss_approx(T[j],n,0,L,U);
            }
            // une fois qu'on a calcule tous les Un[j] pour l'iteration n, on les met dans U[j] pour les reutiliser dans l'iteration suivante
            for (j = 0; j<n; j++) {
                U[j] = Un[j];
            }
        }
    }

    // on libere l'espace memoire des tableaux alloues dynamiquement
    free(T);
    free(U);
    free(Un);

    return u_n;
}