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


/*************************
 fonctions a integrer pour
 les coefficients
**************************/

double B0_f3_f0_noyaux(double x, int m, double alpha);
double Am_f3_noyaux(double x, int m, double alpha);
double Am_f0_noyaux(double x, int m, double alpha);
double Bm_f3_noyaux(double x, int m, double alpha);
double Bm_f0_noyaux(double x, int m, double alpha);
double B0_f3_f0_adomain(double x, int m, double alpha);
double Am_f3_adomain(double x, int m, double alpha);
double Am_f0_adomain(double x, int m, double alpha);
double Bm_f3_adomain(double x, int m, double alpha);
double Bm_f0_adomain(double x, int m, double alpha);

/************************
 fonctions coefficients
*************************/

double A_0_noyaux();
double B_0_noyaux(double alpha);
double A_m_noyaux(int m, double alpha);
double B_m_noyaux(int m, double alpha);
double A_0_adomain();
double B_0_adomain(double alpha);
double A_m_adomain(int m, double alpha);
double B_m_adomain(int m, double alpha);

/************************
 fonction T tilde approchee
*************************/

double T_tilde_noyaux(double x, double y, double alpha);
double T_tilde_adomain(double x, double y, double alpha);






#endif