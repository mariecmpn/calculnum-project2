#ifndef noyaux_h
#define noyaux_h

double k(double x, double t);

double noyaux_iter(double x, double alpha);
double inte_h(double x, int m, double alpha);
double h(double x);
double gauss_approx(double x, int n, double a, double b, double val[5]);

#endif