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
 
int main (int argc, char *argv[]){

   	Random rnd;				//creating a random element to initialize the random generator
	rand_initialize(rnd);



//1. estimation of I with uniform distribution
	cout<<endl<<"//1. estimation of I with uniform distribution"<<endl<<pi-M_PI<<endl;
	ofstream fout1("output2.1.1.txt");
   	int N = 100;				// #blocks
  	int M = 10000;	 			// total throws
  	int L = M/N; 				//throws for each block
  
	double sum=0.;
	double av=0.;
	double av_parz=0.;
	double av2=0.;
	double av2_parz=0.;
	for (int i=0; i<N; i++){
		sum=0.;
		for(int j=0; j<L; j++){
			sum+=(pi/2)*cos((pi/2)*rnd.Rannyu());
		}
		av_parz+=sum/L;
		
		av=av_parz/(i+1);		
		av2_parz+=pow(sum/L,2);
		av2=av2_parz/(i+1);
		fout1<<(i+1)<<"\t"<<fixed<<av<<"\t"<<error(av,av2,i)<<endl;
		cout<<(i+1)<<"\t"<<fixed<<av<<"\t"<<error(av,av2,i)<<endl;
	}
	
	fout1.close();
// 2. estimation of I with importance sampling
	cout<<endl<<"// 2. estimation of I with importance sampling"<<endl;
	ofstream fout2("output2.1.2.txt");
	av=0.;
	av_parz=0.;
	av2=0.;
	av2_parz=0.;
  	for (int i=0; i<N; i++){
		sum=0.;
		for(int j=0; j<L; j++){
			double x=rnd.retta();
			sum+=(pi/2)*cos((pi/2)*x)*(1/(2*(1-x)));
		}
		av_parz+=sum/L;
		
		av=av_parz/(i+1);		
		av2_parz+=pow(sum/L,2);
		av2=av2_parz/(i+1);
		fout2<<(i+1)<<"\t"<<fixed<<av<<"\t"<<error(av,av2,i)<<endl;
		cout<<(i+1)<<"\t"<<fixed<<av<<"\t"<<error(av,av2,i)<<endl;
	}
	
	fout2.close();
  

	



   rnd.SaveSeed();
      
   return 0;
}

