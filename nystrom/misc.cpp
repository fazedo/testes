#include "misc.h"


namespace quadraturas{
	bool verbose=false;

	void Simpson(int N, std::vector<double> *nos, std::vector<double> *pesos, double a, double b){
		if (verbose) std::cout <<"Quadratura: Simpson N = "<<N << std::endl;

		nos->resize(N);
		pesos->resize(N);

		for (int l=0; l<N; l++)	nos->at(l) = a + (b-a) * l / (N-1.0);


		pesos->at(0) = (b - a) / (3.0 * (N-1));
		pesos->at(N-1) = (b - a) / (3.0 * (N-1));
		for (int l=1;l<N-1;l+=2) pesos->at(l) = (b-a) * 4.0/(3.0*(N-1));
		for (int l=2;l<N-2;l+=2) pesos->at(l) = (b-a) * 2.0/(3.0*(N-1));


	}


	void Boole (int N, std::vector<double> *nos, std::vector<double> *pesos, double a, double b) {
		if (N%4!=1) {
			printf("Não existe quadratura de boole com %d pontos!",N);
			return ;
		}
		if (verbose) std::cout <<"Quadratura: Boole N = "<<N << std::endl;

		nos->resize(N);
		pesos->resize(N);


		for (int l=0; l<N; l++)	nos->at(l) = a + (b-a) * l / (N-1.0);


		pesos->at(0) = (b-a)*28/90.0/(N-1.0);

		for (int i=1;i<N;i+=4) {
			pesos->at(i)   = 4.0*16.0*(b-a)/45.0/(N-1.0);
			pesos->at(i+1) = 4.0*12.0*(b-a)/90.0/(N-1.0);
			pesos->at(i+2) = 4.0*16.0*(b-a)/45.0/(N-1.0);
			pesos->at(i+3) = 4.0*7.0*(b-a)/45.0/(N-1.0);
		}
		pesos->at(N-1) = 28*(b-a)/90.0/(N-1.0);

	}



void GaussLegendre(int N, std::vector<double> *nos, std::vector<double> *pesos, double a, double b){
		double node, weight;

		if (verbose) std::cout <<"Quadratura: Gauss-Legendre N = "<< N << std::endl;

		nos->resize(N);
		pesos->resize(N);

		gsl_integration_glfixed_table * table = gsl_integration_glfixed_table_alloc (N);

		for (int i=0; i<N;i++) {
			gsl_integration_glfixed_point (a, b, i, &node, &weight, table);
			pesos->at(i) = weight;
			nos->at(i)   = node;

		}
		return;
		}

}
