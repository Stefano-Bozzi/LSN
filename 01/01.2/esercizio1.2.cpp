/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>

using namespace std;
 
int main (int argc, char *argv[]){
   	Random rnd;				//creating a random element to initialize the random generator
   	rand_initialize(rnd);
		
   	vector <int> N{1,2,10,100};
   	int n=10000;
   	vector <double> S{0., 0., 0.};
   	vector <string> filename (3);
   	vector <string> type {"std","exp","Cauchy"};
   	
   	for (int i=0; i<4; i++){
   		for(int j=0; j<3; j++) filename[j] = type [j]+to_string(N[i]) + ".txt";
   		ofstream fouts(filename[0]);
   		ofstream foute(filename[1]);
   		ofstream foutc(filename[2]);
   		
   		for (int k=0; k<n; k++){
   			S={0.,0.,0.};
   			for (int l=0; l<N[i];l++){
   				S[0]+=rnd.Rannyu();
   				S[1]+=rnd.Exp(1.);
   				S[2]+=rnd.Cauchy(1.,0.);
   			}
   			fouts<<S[0]/N[i]<<endl;
   			foute<<S[1]/N[i]<<endl;
   			foutc<<S[2]/N[i]<<endl;
   		}
   		fouts.close();
   		foute.close();
   		foutc.close();
   	}
   	
	
   	rnd.SaveSeed();
   	return 0;
}	

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
