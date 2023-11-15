#ifndef fonctions_h
#define fonctions_h


/************************
 fonctions conditions 
 aux limites et exactes
*************************/

double T_ex(double x, double y);
double f_0(double x, int m, double alpha);
double f_3_ex(double x);
double q_0(double x);

/************************
 fonctions a approcher
*************************/

//double inte_h(double x, int m, double alpha);
//double h(double x);
double inte_f3(double x, int m, double alpha);
double f_3(double x, double alpha);

/*************************
 fonctions a integrer pour
 les coefficients
**************************/

double B0_f3_f0(double x, int m, double alpha);
double Am_f3(double x, int m, double alpha);
double Am_f0(double x, int m, double alpha);
double Bm_f3(double x, int m, double alpha);
double Bm_f0(double x, int m, double alpha);

/************************
 fonctions coefficients
*************************/

double A_0();
double B_0(double alpha);
double A_m(int m, double alpha);
double B_m(int m, double alpha);

/************************
 fonction T tilde approchee
*************************/

double T_tilde(double x, double y, double alpha);






#endif