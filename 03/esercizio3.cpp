#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;
using std::fixed;

double pi=acos(-1);
double error(double av, double av2, int n){	//function for calculate the error
	if(n==0)
		return 0;
	else
		return sqrt((av2- av*av)/n);
	}
 double max(double a,double b){
 		if (a>=b) return a;
 		else return b;
 	}
int main (int argc, char *argv[]){

   	Random rnd;				//creating a random element to initialize the random generator
	rand_initialize(rnd);
	int M=10E4; 		//tot di simulazioni
	int N=100;		//numero di cicli
	int L=M/N;		//numero di simulazioni per ciclo
	int n = 100;		//numero di intervalli in cui suddividere il tempo T-t_0
	double t_0=0.;
	double T=1.;
	double dt=(T-t_0)/n;
	double S_T=0.;
	double S_0=100.;
	double K=100.;
	double r=0.1;
	double vol=0.25;
	double sum_C=0.;
	double sum_P=0.;
	double C=0.;
	double P=0.;
	double C_2=0.;
	double P_2=0.;
	//vector <double> C(N);
	//vector <double> P(N);
	
	/*for (int i=0; i<N; i++){
		sum_C=0.;
		sum_P=0.;
	}*/
	
	ofstream fout1C("output1C.txt");
	ofstream fout1P("output1P.txt");
	for (int i=0;i<N;i++){
		sum_C=0.;
		sum_P=0.;
		for (int j=0; j<L; j++){
			S_T=S_0*exp((r-((pow(vol,2))/2))*T+vol*rnd.Gauss(0,1)*sqrt(T));
			sum_C+=exp(-r*T)*max(0,S_T-K);
			sum_P+=exp(-r*T)*max(0,K-S_T);
		}
		C+=sum_C/L;
		P+=sum_P/L;
		C_2+=pow(sum_C/L,2);
		P_2+=pow(sum_P/L,2);
		fout1C<<(i+1)<<"\t"<<C/(i+1)<<"\t"<<error(C/(i+1),C_2/(i+1),i)<<endl;
		fout1P<<(i+1)<<"\t"<<P/(i+1)<<"\t"<<error(P/(i+1),P_2/(i+1),i)<<endl;
		
		
	}
	fout1C.close();
	fout1P.close();
	
//DISCRETIZZATO

	ofstream fout2C("output2C.txt");
	ofstream fout2P("output2P.txt");
	C=0.;
	P=0.;
	C_2=0.;
	P_2=0.;
	for (int i=0;i<N;i++){
		sum_C=0.;
		sum_P=0.;
		for (int j=0; j<L; j++){
			S_T=S_0;
			for (int t=0; t<n; t++){
				S_T=S_T*exp((r-((pow(vol,2))/2))*dt+vol*rnd.Gauss(0,1)*sqrt(dt));
			}
			sum_C+=exp(-r*T)*max(0,S_T-K);
			sum_P+=exp(-r*T)*max(0,K-S_T);
		}
		C+=sum_C/L;
		P+=sum_P/L;
		C_2+=pow(sum_C/L,2);
		P_2+=pow(sum_P/L,2);
		fout2C<<(i+1)<<"\t"<<C/(i+1)<<"\t"<<error(C/(i+1),C_2/(i+1),i)<<endl;
		cout<<(i+1)<<"\t"<<C/(i+1)<<"\t"<<error(C/(i+1),C_2/(i+1),i)<<endl;
		fout2P<<(i+1)<<"\t"<<P/(i+1)<<"\t"<<error(P/(i+1),P_2/(i+1),i)<<endl;
		
		
	}
	fout2C.close();
	fout2P.close();	
	rnd.SaveSeed();
      
   	return 0;
}
	
