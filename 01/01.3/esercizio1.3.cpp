#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;
using std::fixed;

double error(double av, double av2, int n){	//function for calculate the error
	if(n==0)
		return 0;
	else
		return sqrt((av2- av*av)/n);
	}
 
int main (int argc, char *argv[]){

   	Random rnd;				//creating a random element to initialize the random generator
	rand_initialize(rnd);



//1. estimation of <pi>
	cout<<endl<<"//1. estimation of <pi>"<<endl;
	ofstream fout1("pi.txt");
   	int N = 100;				// #blocks
  	int M = 10000000;			// total throws          
  	int N_thr = M/N; 				//throws for each block
  	int N_hit = 0;
  	double L=0.3;
  	double d=1.;
  	
	double pigreco=0.;
	double pigreco_parz=0.;
	double pigreco2=0.;
	double pigreco2_parz=0.;
	double y=0;
	double theta=0;
	
	for (int i=0; i<N; i++){
		N_hit=0;
		for(int j=0; j<N_thr; j++){
			y= rnd.Rannyu(0,d);
			theta= rnd.Angle();
			if(y+L*sin(theta)>=d || y+L*sin(theta)<=0.0){N_hit+=1;}
			
		}
		pigreco_parz+=2*L*N_thr/(N_hit*d);
		pigreco2_parz+=pow(2*L*N_thr/(N_hit*d),2);
		
		pigreco=pigreco_parz/(i+1);		
		pigreco2=pigreco2_parz/(i+1);
		fout1<<(i+1)<<"\t"<<fixed<<pigreco<<"\t"<<error(pigreco,pigreco2,i)<<endl;
		cout<<(i+1)<<"\t"<<fixed<<pigreco<<"\t"<<error(pigreco,pigreco2,i)<<endl;
	}
	
	fout1.close();

   rnd.SaveSeed();
      
   return 0;
}

