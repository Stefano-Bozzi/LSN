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



//1. estimation of <r>
	cout<<endl<<"//1. estimation of <r>"<<endl;
	ofstream fout1("output1.01.1.txt");
   	int N = 100;				// #blocks
  	int M = 100000; 			// total throws
  	int L = M/N; 				//throws for each block
  
	double sum=0.;
	double av=0.;
	double av_parz=0.;
	double av2=0.;
	double av2_parz=0.;
	for (int i=0; i<N; i++){
		sum=0.;
		for(int j=0; j<L; j++){
			sum+=rnd.Rannyu();
		}
		av_parz+=sum/L;
		
		av=av_parz/(i+1);		
		av2_parz+=pow(sum/L,2);
		av2=av2_parz/(i+1);
		fout1<<(i+1)<<"\t"<<fixed<<av<<"\t"<<error(av,av2,i)<<endl;
		cout<<(i+1)<<"\t"<<fixed<<av<<"\t"<<error(av,av2,i)<<endl;
	}
	
	fout1.close();
  	

//2. estimation of var
	cout<<endl<<"//2. estimation of var"<<endl;
	ofstream fout2("output1.01.2.txt");
	av=0.;
	av_parz=0.;
	av2=0.;
	av2_parz=0.;
  	for (int i=0; i<N; i++){
		sum=0.;
		for(int j=0; j<L; j++){
			sum+=pow((rnd.Rannyu()-0.5),2);		//nerly the same process used for estimation of <r>
		}
		av_parz+=sum/L;
		
		av=av_parz/(i+1);		
		av2_parz+=pow(sum/L,2);
		av2=av2_parz/(i+1);
		fout2<<(i+1)<<"\t"<<fixed<<av<<"\t"<<error(av,av2,i)<<endl;
		cout<<(i+1)<<"\t"<<fixed<<av<<"\t"<<error(av,av2,i)<<endl;
	}
	
	fout2.close();
  

//3. chi2 test
	cout<<endl<<"//3. chi2 test"<<endl;
	ofstream fout3("output1.01.3.txt");
	M=100;
	N=10000;
	vector <int> o;
	for (int i=0; i<M; i++){
		o.push_back(0);				//create vector for counting probability
	}
	for (int l=0; l<100; l++){			//repeat the calculation of chi2 100 times
		for (int i=0; i<M; i++) o[i]=0;		
		double chi2=0.;
		for (int i=0; i<N; i++){
			int classNumber=rnd.Rannyu()*M;
			++o[classNumber];		//counting probability
		}
		for (int j=0; j<M; j++){
			chi2+=pow((o[j]-(N/M)),2)/(N/M);//calculating chi2
		}
		fout3<<chi2<<endl;
		cout<<chi2<<endl;
	}
		fout3.close();

   rnd.SaveSeed();
      
   return 0;
}

