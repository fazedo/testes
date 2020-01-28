#include <iostream>
#include <vector>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "misc.h"

double kernel (double x, double y);
double fonte (double x);
double nystron (int N, std::vector<double> nos, std::vector<double> pesos);
double interpola(int N, double sigma, std::vector<double> pesos, std::vector<double> nos, gsl_vector* u, double x);
double Kint (double x);

int main(){
  int N = 17; //numero de pontos
  double a = 0;
  double b = 1;

  std::vector<double> nos, pesos;


  //quadraturas::Simpson(N, &nos, &pesos, a, b) ;
  //quadraturas::Boole(N, &nos, &pesos, a, b) ;
  quadraturas::GaussLegendre(N, &nos, &pesos, a, b) ;



  double sigma = .7;


  gsl_matrix* matriz; // matrix principal
  matriz = gsl_matrix_alloc(N, N);


  for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if ( i != j ){
    gsl_matrix_set(matriz, i, j,  - sigma * pesos.at(j) * kernel (nos.at(i), nos.at(j)));
  }

  for (int i = 0; i < N; i++){
    double x = 1;
    for (int l = 0; l < N; l++) if (l != i) x += sigma * pesos.at(l) * kernel(nos.at(l), nos.at(i));
    x -= sigma * Kint(nos.at(i));

    gsl_matrix_set(matriz, i, i,  x);
  }


  //constroi vetor
  gsl_vector* f;
  f = gsl_vector_alloc(N);
  for (int i = 0; i < N; i++) gsl_vector_set(f, i, fonte(nos.at(i)));

  //resolve sistema
  gsl_vector* u;
  u = gsl_vector_alloc(N);
  int s;
  gsl_permutation * p = gsl_permutation_alloc (N);
  gsl_linalg_LU_decomp (matriz, p, &s);
  gsl_linalg_LU_solve (matriz, p, f, u);

  for (int i = 0; i<N; i++) std::cout << "x = " << nos.at(i) << " \t u(x) = " <<  gsl_vector_get(u, i) << "  " << interpola(N, sigma, pesos, nos, u, nos.at(i)) << std::endl;


  return 0;
}



double kernel(double x, double y){

	return -log(fabs(x-y));

}

//int( k(x,y), y=Left..Right)
double Kint (double x){
	if ((x<=0) || (x>=1)) return 1;

	return -(-1.0 + x*log(x) + (1-x)*log(1-x)) ;


}


double fonte (double x) {

  return 1.0 -  Kint(x);


}



double interpola(int N, double sigma, std::vector<double> pesos, std::vector<double> nos, gsl_vector* u, double x){
  double xx = 0;
  double yy = 0;

 for (int j = 0; j < N; j++) {
   double xi = nos.at(j);
   if (x != xi) {
     xx += pesos.at(j) * kernel (x, nos.at(j)) * gsl_vector_get(u,j);
     yy += pesos.at(j) * kernel (x, nos.at(j));
   }
  }

  return (sigma * xx + fonte(x)) / (1 + sigma * yy - sigma * Kint(x));
}
