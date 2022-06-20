#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>
#include <cstdlib>
#include "lib.h"

using namespace std;

const int x=0;
const int y=1;
const int L2=0;


void CheckCoincidence(vector<vector<double>> &v, int N, Random &rnd){
	int flag=0;
	do{
		flag=0;
		for (int i=1; i<N+1; i++){			//un ciclo for per vedere se tutte le città sono diverse tra loro
			for (int j=i+1; j<N+1; j++){
				if (v[j][x]==v[i][x]&&v[j][y]==v[i][y]){
					cout<<"AAAAAAAAAAAAAAAAAAAAAAAAA"<<endl;
					flag=1;
					v[j][x]=rnd.Rannyu();
					v[j][y]=rnd.Rannyu();
				}
			}
		}
	}while(flag==1);	
}


void WriteBestPath(vector<vector<vector<double>>> Pop, int N, ofstream &cityPath1, const std::string& name, int n){
	for (int j=1; j<N+1; j++){
					//cout<<"writing"<<name<<"_"<<to_string(n)<<".dat"<<endl;
					cityPath1.open(name+"_"+to_string(n)+".dat",ios::app);
					cityPath1<<j<<"\t"<<Pop[0][j][x]<<"\t"<<Pop[0][j][y]<<"\n";
					cityPath1.close();
	}
	cityPath1.open(name+"_"+to_string(n)+".dat",ios::app);
					cityPath1<<1<<"\t"<<Pop[0][1][x]<<"\t"<<Pop[0][1][y]<<"\n";
					cityPath1.close();
}


void CalculateL2(vector<vector<double>> &v, int N){
	double Dij;
	vector<double> result;
		
	for(int i=1; i<N; i++){
		Dij = pow((v[i][x] - v[i+1][x]), 2) + pow((v[i][y] - v[i+1][y]), 2);
		result.push_back(Dij);
	}
	Dij = pow((v[N][x] - v[1][x]), 2) + pow((v[N][y] - v[1][y]), 2);
	result.push_back(Dij);
	v[L2][0]=0.;
	for(int i=1; i<N+1; i++){
		v[L2][0]=v[L2][0]+result[i];
	}
	//cout<<"ho calcolato L2="<<v[L2][0]<<endl;
}


void mescola(vector<vector<double>> &v, int N, Random &rnd){
	int pos;
	for (int i=2; i<N+1; i++){
		pos=(int)(rnd.Rannyu(2,N+1));
		v[i].swap(v[pos]);
	}
	CalculateL2(v,N);
	/*cout<<"ho MESCOLATO=:"<<endl;
	for (int j=1; j<N+1; j++){
			cout<<"|"<<v[j][x]<<"\t"<<";"<<"\t"<<v[j][y]<<"|";
			cout<<endl;
	}*/
}


void ordina(vector<vector<vector<double>>> &v, int N, int Ndim){	//funzione per ordinare la Popolazione secondo fitness decrescente(if vi>vj) se voglio crescente(if vi<vj)
	for (int i=0; i<Ndim; i++){					//nella selezione dovrò selezionare gli individui che sono ai posti più bassi [0,1,2,...]
		for (int j=i+1; j<Ndim; j++){
			if (v[i][L2][0]>v[j][L2][0]){
				v[i].swap(v[j]);
			}
		}
	}
	/*cout<<"ho ORDINATO=:"<<endl;
	for (int i=0; i<Ndim; i++){
		cout<<"\n\n\n\n"<<"POPOLAZIONE NUMERO 1, lista di città n:"<<i+1<<"\n"<<"somma iniziale= "<<v[i][0][0]<<endl;
		cout<<"| 	x 	; 	y	 |"<<endl;
		for (int j=1; j<N+1; j++){
				cout<<"|"<<v[i][j][x]<<"\t"<<";"<<"\t"<<v[i][j][y]<<"|";
				cout<<endl;
		}
	}*/
}


int seleziona(int Ndim, Random &rnd){
	int j;
	double r = rnd.Rannyu(0., 1.);
	j = Ndim*pow(r, 3);
return j;
}


void permCoppie(vector<vector<double>> &v, int N, Random &rnd){
	int a,b;
	a=(int)(rnd.Rannyu(2,N+1));
	b=(int)(rnd.Rannyu(2,N+1));
	v[a].swap(v[b]);
	CalculateL2(v,N);
	//cout<<a<<"\t"<<b<<endl;
	/*cout<<"ho PERMUTATO A COPPIE=:"<<endl;
	for (int j=1; j<N+1; j++){
			cout<<"|"<<v[j][x]<<"\t"<<";"<<"\t"<<v[j][y]<<"|";
			cout<<endl;
	}*/
}


void permContig(vector<vector<double>> &v, int N, Random &rnd){
	int a,b,m;
	m=(int)(rnd.Rannyu(1,N/2));//numero di celle da cambiare
	a=(int)(rnd.Rannyu(2,N+2-2*m));//punto di inizio per celle nella prima metà
	b=(int)(rnd.Rannyu(a+m,N+2-m));//punto di inizio per celle nella seconda metà
	for(int i=a; i<a+m; i++){
		v[i].swap(v[b]);
		b++;
	}
	CalculateL2(v,N);
	/*cout<<"ho PERMUTATO CONTIGUAMENTE=:"<<endl;
	for (int j=1; j<N+1; j++){
			cout<<"|"<<v[j][x]<<"\t"<<";"<<"\t"<<v[j][y]<<"|";
			cout<<endl;
	}*/
}


void inversione(vector<vector<double>> &v, int N, Random &rnd){
	int a,m;
	m=(int)(rnd.Rannyu(2,N));//numero di celle da invertire
	a=(int)(rnd.Rannyu(2,N+2-m));//punto di inizio di inversione di celle
	for (int i=0; i<m/2; i++){
		v[a+i].swap(v[a+m-1-i]);
	}
	CalculateL2(v,N);
	/*cout<<"ho INVERTITO=:"<<endl;
	for (int j=1; j<N+1; j++){
			cout<<"|"<<v[j][x]<<"\t"<<";"<<"\t"<<v[j][y]<<"|";
			cout<<endl;
	}*/
}


void shiftContig(vector<vector<double>> &v, int N, Random &rnd){
	int m,n,start;
	m=(int)(rnd.Rannyu(1,N-1));
	n=(int)(rnd.Rannyu(1,N-m));	
	start=(int)(rnd.Rannyu(2,2+N-m-n));
	vector<vector<double>> copia;
	copia=v;
	for (int i=start; i<start+m; i++){
		v[i+n] = copia[i];
	}
	for (int i=start; i<start+n; i++){
		v[i] = copia[i+m];
	}
	CalculateL2(v,N);
	/*cout<<"ho SCALATO=:"<<endl;
	for (int j=1; j<N+1; j++){
			cout<<"|"<<v[j][x]<<"\t"<<";"<<"\t"<<v[j][y]<<"|";
			cout<<endl;
	}*/
}
	

void crossover(vector<vector<double>> &a, vector<vector<double>> &b, int N, Random &rnd){
 	int m,k;
	m=(int)(rnd.Rannyu(2,N+1));
	k=m;
	//cout<<k<<endl;
	vector<vector<double>> a_copia, b_copia;
	a_copia=a;
	b_copia=b;
	for (int i=2; i<N+1; i++){		//  ciclo per il cabiamento di a con ordine di b
		for (int j=m; j<N+1; j++){
			if(b[i]==a_copia[j]){	
				a[k]=b[i];	
				k++;		
			}			
		}				
	}
	k=m;
	for (int i=2; i<N+1; i++){		//  ciclo per il cabiamento di a con ordine di b
		for (int j=m; j<N+1; j++){
			if(a_copia[i]==b_copia[j]){	
				b[k]=a_copia[i];	
				k++;		
			}			
		}				
	}
	CalculateL2(a,N);
	CalculateL2(b,N);
	/*cout<<"ho FATTO IL CROSSOVER=:"<<endl;
	for (int j=1; j<N+1; j++){
			cout<<"|"<<a[j][x]<<"\t"<<";"<<"\t"<<a[j][y]<<"|";
			cout<<endl;
	}cout<<endl;
	for (int j=1; j<N+1; j++){
			cout<<"|"<<b[j][x]<<"\t"<<";"<<"\t"<<b[j][y]<<"|";
			cout<<endl;
	}*/
}


vector<vector<vector<double>>> GeneraPopStart(vector<vector<double>> v, int N, int Ndim, Random &rnd){
	CalculateL2(v, N);
	vector<vector<vector<double>>> Pop;
	Pop.push_back(v);
	for (int i=0; i<Ndim ;i++){
		mescola(v, N, rnd);
		Pop.push_back(v);
		CalculateL2(Pop[i], N);
	}
	return Pop;
}


void NewPop(vector<vector<vector<double>>> &Pop_old, int N, int Ndim, double p_c, double p_m, Random &rnd){
	int gen1, gen2;
	ordina(Pop_old, N, Ndim);
	vector<vector<vector<double>>> Pop_New, Pop_copia;
	Pop_copia=Pop_old;
	for(int l=0; l<Ndim/2; l++){
		do{
			gen1=seleziona(Ndim, rnd);
			gen2=seleziona(Ndim, rnd);
		}while(gen1 == gen2);
				
		if (rnd.Rannyu()<p_c){
			crossover(Pop_old[gen1], Pop_old[gen2], N, rnd);
		}
		
		if (rnd.Rannyu()<p_m){
			permCoppie(Pop_old[gen1], N, rnd);
		}
		
		if (rnd.Rannyu()<p_m){
			permCoppie(Pop_old[gen2], N, rnd);
		}
		if (rnd.Rannyu()<p_m){
			permContig(Pop_old[gen1], N, rnd);
		}
		
		if (rnd.Rannyu()<p_m){
			permContig(Pop_old[gen2], N, rnd);
		}
		if (rnd.Rannyu()<p_m){
			inversione(Pop_old[gen1], N, rnd);
		}
		
		if (rnd.Rannyu()<p_m){
			inversione(Pop_old[gen2], N, rnd);
		}
		if (rnd.Rannyu()<p_m){
			shiftContig(Pop_old[gen1], N, rnd);
		}
		
		if (rnd.Rannyu()<p_m){
			shiftContig(Pop_old[gen2], N, rnd);
		}
		Pop_New.push_back(Pop_old[gen1]);
		Pop_New.push_back(Pop_old[gen2]);
		Pop_old=Pop_copia;
		//cout<<Pop_New.size()<<endl;
	}
	Pop_old=Pop_New;
	for (int i=0; i<Ndim ;i++){
		CalculateL2(Pop_old[i], N);		
	}
}
