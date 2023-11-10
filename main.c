#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "donnees.h" 
#include "methodesnum.h"


int main() {

    /***************************
     definition des variables 
    ****************************/ 

    int i, j; // entiers pour les boucles for
    int Nx = 51; // nombre de points de maillage en x
    int Ny = 51; // nombre de points de maillage en y
    int Napp = Nx, Nex = Nx;
    double pasx = 1./(Nx-1); // pas de maillage en x
    double pasy = 1./(Ny-1); // pas de maillage en x
    double *Points = malloc(sizeof(int[Nx])); // tableau qui contient les valeurs des x_i
    double *Y = malloc(sizeof(int[Ny])); // tableau qui contient les y_i

    double Alpha[22]; // tableau qui contient les differentes valeurs de alpha

    // recuperation des donnees du probleme definies dans donnees.c
    double L = recup_L(L), H = recup_H(H);
    int M = recup_M(M), n = recup_n(n);

    double *T_app = malloc(sizeof(int[Napp])); // tableau qui contient les images des x_i par T_alpha calculees numeriquement
    double *T_exact = malloc(sizeof(int[Nex])); // tableau qui contient les images des x_i par T_ex calculees numeriquement

    double y = H; // on se place sur Gamma_3



    /***************************
     remplissage des tableaux
    ****************************/ 

    // remplissage de Points
    /*double x_i = 0.;
    for (i = 0; i<N; i++) {
        Points[i] = x_i;
        x_i  = x_i+pas;
        //printf("%f\n", Points[i]);
    }*/

    // remplissage de Alpha
    Alpha[0] = 10.;
    Alpha[1] = 1.;
    for (i = 2; i<22; i++) {
        Alpha[i] = 1./pow(10.,i+1);
    }


    // remplissage de Points
    double x_i = 0.;
    //double y_i = 0.;
    for (i = 0; i<Nx; i++) {
        //Y[i] = y_i;
        //y_i  = y_i+pas;
        Points[i] = x_i;
        x_i  = x_i+pasx;
        //printf("%f", Points[i]);
    }


    /***************************
     Calcul pour plusieurs alpha
     sur Gamma_3 de la solution
     approchee
    ****************************/ 

    y = H; // on se place sur Gamma_3
    FILE *approche;
    approche = fopen("approche_1.txt", "w"); // on ouvre le fichier en ecriture
    for (i = 0; i < 22; i++) {
        for (j = 0; j < Nx; j++) {
            T_app[j] = f_3(Points[j], Alpha[i]); // on calcule T_tilde pour alpha et x_i: sur Gamma_3 T_tilde = f_3 Alpha[i]
            //printf("%f", Points[j]);
            fprintf(approche, "%g", T_app[j]); // on l'enregistre dans le fichier approche_1.txt
            fputs(" ", approche);
        }
        fputs("\n", approche); // on change de ligne quand on change de alpha
    }

    fclose(approche); // on ferme le fichier une fois qu'il est rempli

    /***************************
    Calcul de la solution exacte
    sur Gamma_3
    ****************************/ 

    FILE *exact;
    exact = fopen("exact_1.txt", "w"); // on ouvre le fichier
    for (i = 0; i<Nx; i++) {
        T_exact[i] = f_3_ex(Points[i]); // calcul de T sur Gamma_3 
        //T_exact[i] = cosh(M_PI*H)*cos(M_PI*Points[i]);
        //printf("%f%s", Points[i], " ");
        //printf("%f\n", T_exact[i]);
        fprintf(exact, "%g", T_exact[i]); // on ecrit dans le fichier
        fputs(" ", exact);
    }

    fclose(exact); // on ferme le fichier

    /***************************
    Calcul des solutions sur le
    domaine Omega pour alpha
    optimal
    ****************************/ 

    double alpha_optim = 1.;
    FILE *app_omega;
    app_omega = fopen("solapp_omega.txt", "w");
    FILE *ex_omega;
    ex_omega = fopen("solex_omega.txt", "w");
    FILE *erreur;
    erreur = fopen("erreur_omega.txt", "w");

    // on calcule les solutions exactes et approchees pour tous les points de maillage du domaine
    double y_i = 0.;
    x_i = 0.;
    for (j = 0; j < Ny; j++) {
        Y[j] = y_i;
        //printf("%f\n", Y[j]);
        x_i = 0.;
        for (i = 0; i < Nx; i++) {
            Points[i] = x_i;
            T_app[i] = T_tilde(x_i, y_i, alpha_optim); 
            T_exact[i] = T_ex(x_i, y_i);
            //T_exact[i] = cosh(M_PI*y_i)*cos(M_PI*x_i);
            fprintf(app_omega, "%g", T_app[i]);
            fputs(" ", app_omega);
            fprintf(ex_omega, "%g", T_exact[i]);
            fputs(" ", ex_omega);
            fprintf(erreur, "%lf", fabs(T_exact[i]-T_app[i]));
            fputs(" ", erreur);
            //printf("%f\n", Points[i]);
            x_i  = x_i+pasx;
        }
        fputs("\n", app_omega); // on change de ligne quand on change de y_i
        fputs("\n", ex_omega);
        fputs("\n", erreur);
        y_i  = y_i+pasy;
    }


    // on ferme les fichiers
    fclose(app_omega); 
    fclose(ex_omega);
    fclose(erreur);


    // on desalloue l'espace memoire des tableaux alloues dynamiquement
    free(Points);
    free(Y);
    free(T_exact);
    free(T_app);

    // on retourne 0
    return 0;
}