#ifndef MISC
#define MISC


#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <vector>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <algorithm>


namespace quadraturas{
	void Simpson(int N, std::vector<double> *nos, std::vector<double> *pesos, double a = -1, double b = 1);
	void GaussLegendre(int N, std::vector<double> *nos, std::vector<double> *pesos, double a = -1, double b = 1);
	void Boole(int N, std::vector<double> *nos, std::vector<double> *pesos, double a = -1, double b = 1);
	double Error (int N, std::vector<double> *v1, std::vector<double> *v2);
}

#endif
