#ifndef adomain_h
#define adomain_h

double k(double x, double t);

double h(double x);

double inte_h(double x, int m, double alpha);

double gauss_approx(double x, int n, double a, double b, double val[5]);

double Adomain(double x, double alpha);

#endif
