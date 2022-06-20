#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>
#include <cstdlib>
#include "Posizione.h"

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
	int N = 100;				//numero di step in un RW
	int N_t = 10000;			//numero di RW totali
	int N_b = 100;				//numero di blocchi
	int L = N_t/N_b;			//numero di RW per blocco
	vector <double> sum(N);
	vector <double> av(N);
	vector <double> av2(N);
	Posizione RW_lat;
	Posizione RW_cont;
	
//DISCRETO
	ofstream fout1("output2.2.discr.txt");
	for (int i=0; i<N; i++) {
		av[i]=0;
		av2[i]=0;
		sum[i]=0;
	}
	for (int j=0; j<N_b; j++){			//ciclo sui blocchi
		for (int i=0; i<N; i++) sum[i]=0;	
		
		for (int i=0; i<L; i++) {		//ciclo sui RW di un blocco
			RW_lat.setX(0.);
			RW_lat.setY(0.);
			RW_lat.setZ(0.);
			for (int k=0; k<N; k++) {	//ciclo sugli step di un RW
				double r=rnd.Rannyu();
				if(r>5./6.){
					RW_lat.setX(RW_lat.getX()+1.);
				}else if(r>4./6.&&r<5./6.){
					RW_lat.setX(RW_lat.getX()-1.);
				}else if(r>3./6.&&r<4./6.){
					RW_lat.setY(RW_lat.getY()+1.);
				}else if(r>2./6.&&r<3./6.){
					RW_lat.setY(RW_lat.getY()-1.);
				}else if(r>1./6.&&r<2./6.){
					RW_lat.setZ(RW_lat.getZ()+1.);
				}else{
					RW_lat.setZ(RW_lat.getZ()-1.);
				}
				sum[k]+= pow(RW_lat.getR(),2);
			}//qui ho ottenuto un RW completo
		}//qui ho fatto un tot di RW e ho riempito un blocco
		
		for (int k=0; k<N; k++){			
			av[k]+=sqrt(sum[k]/L);
			av2[k]+=sum[k]/L;
		}//qui ho trovato il valore radq(media(mod(r)^2)) e il suo quadrato e li sommo così da ottenere poi un tot di quantità su cui fare la media totale e la varianza totale
	}
	for (int k=0; k<N; k++){			
			fout1<<k+1<<"\t"<<fixed<<av[k]/N_b<<"\t"<<error(av[k]/N_b,av2[k]/N_b,N_b)<<endl;
			cout<<k+1<<"\t"<<fixed<<av[k]/N_b<<"\t"<<error(av[k]/N_b,av2[k]/N_b,N_b)<<endl;
		}
	fout1.close();
	
//CONTINUO
	ofstream fout2("output2.2.cont.txt");
	for (int i=0; i<N; i++) {
		av[i]=0;
		av2[i]=0;
		sum[i]=0;
	}
	for (int j=0; j<N_b; j++){			//ciclo sui blocchi
		for (int i=0; i<N; i++) sum[i]=0;	
		
		for (int i=0; i<L; i++) {		//ciclo sui RW di un blocco
			RW_cont.setX(0.);
			RW_cont.setY(0.);
			RW_cont.setZ(0.);
			for (int k=0; k<N; k++) {	//ciclo sugli step di un RW
				double r=rnd.Rannyu();
				double theta=acos(1-2*r);
				double phi=rnd.Rannyu(0,2*pi);
				
				RW_cont.setX(RW_cont.getX()+sin(theta)*cos(phi));
				RW_cont.setY(RW_cont.getY()+sin(theta)*sin(phi));
				RW_cont.setZ(RW_cont.getZ()+cos(theta));
				
				sum[k]+= pow(RW_cont.getR(),2);
			}//qui ho ottenuto un RW completo
		}//qui ho fatto un tot di RW e ho riempito un blocco
		
		for (int k=0; k<N; k++){			
			av[k]+=sqrt(sum[k]/L);
			av2[k]+=sum[k]/L;
		}//qui ho trovato il valore radq(media(mod(r)^2)) e il suo quadrato e li sommo così da ottenere poi un tot di quantità su cui fare la media totale e la varianza totale
	}
	for (int k=0; k<N; k++){		
			fout2<<k+1<<"\t"<<fixed<<av[k]/N_b<<"\t"<<error(av[k]/N_b,av2[k]/N_b,N_b)<<endl;
			cout<<k+1<<"\t"<<fixed<<av[k]/N_b<<"\t"<<error(av[k]/N_b,av2[k]/N_b,N_b)<<endl;
		}
	
	fout2.close();	
	
	
   rnd.SaveSeed();
      
   return 0;
}

