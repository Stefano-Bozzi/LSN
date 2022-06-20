//HO FATTO UN CODICE CHE PRENDE LE CITTÀ IN ENTRATA E CI FA UNA LISTA UGUALE PER TUTTI I RANK E MESCOLA LORDINE DEL VETTORE RANK UGUALE PER TUTTI

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>
#include <cstdlib>
#include "lib.h"
#include "mpi.h"
using namespace std;

int x=0;
int y=1;
int L2=0;

	
int main(int argc, char* argv[]){

	int size, rank;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	int N=0;
	int Ngen=1000;
	int Ndim=200;
	int Nmigr=5;
	double p_c=50./100.;
	double p_m=10/100.;
	double avL2=0;
	string namefile;

	ofstream citypathfinal, citypathstart, averegeL2, bestaveregeL2;

	Random rnd;				//random element unico per ogni processo (evoluzione dei sistemi separata e indipendente)
	Random rnd_com;				//random element comune a tutti i processi, serve a generare liste per scambio tra processi 
	
//Random rnd/rnd_com initialization
	int seed0[4];
	int p01, p02;
	ifstream Primes0("Primes");
	if (Primes0.is_open()){
		Primes0 >> p01 >> p02;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes0.close();
	
	ifstream input0("seed.in");
	string property0;
	if (input0.is_open()){
		while ( !input0.eof() ){
			input0 >> property0;
			if( property0 == "RANDOMSEED" ){
		         	input0 >> seed0[0] >> seed0[1] >> seed0[2] >> seed0[3];
		         	rnd_com.SetRandom(seed0,p01,p02);
		      	}
		}
		input0.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	for(int i=0; i<size; i++){
		if(rank == i){
			ifstream input("seed_" + to_string(i) + ".in");
			string property;
			if (input.is_open()){
	   			while ( !input.eof() ){
	      			input >> property;
	      			if( property == "RANDOMSEED" ){
	        			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	        			rnd.SetRandom(seed,p1,p2);
	      			}
	   			}
	  			input.close();
			} else cerr << "PROBLEM: Unable to open seed.in" << endl;
			
		}
	}
//fine inizializzazione random

//lettura da file
	ifstream ifile;
	ifile.open("American_capitals.in");
	vector<double> city{0.,0.};
	vector<double>initialL2{0.};
	vector<vector<double>> capitals;
	capitals.push_back(initialL2);
	//check to see that the file was opened correctly:
	if (!ifile.is_open()) {
		cerr << "There was a problem opening the input file: American_capitals.in!\n";
		exit(1);
	}
	double longitude = 0.0;
	double latitude = 0.0;
	string word;
	
	ifile >> word>> word>> word>> word;
	while (ifile >> word >> word >> longitude >> latitude) {
		city[x]=longitude;
		city[y]=latitude;
		capitals.push_back(city);
	}

	N=capitals.size()-1;

	MPI_Barrier(MPI_COMM_WORLD);
	vector<vector<vector<double>>> Pop;
	Pop = GeneraPopStart(capitals, N, Ndim, rnd);
	ordina(Pop, N, Ndim);
	namefile="start_Path_";
	namefile=namefile+to_string(size);
	WriteBestPath(Pop, N, citypathfinal, namefile ,rank);
	vector<int> v_rank;
	
	if (rank==0)cout<<rank<<"\tv_rank:\n";
	for (int i=0; i<size; i++){
		v_rank.push_back(i);
		if (rank==0)cout<<rank<<"\t"<<v_rank[i]<<endl;
	}
	for (int i=0; i<size; i++){
		int j=(int)(rnd_com.Rannyu(0,size));
		int copy=v_rank[i];
		v_rank[i]=v_rank[j];
		v_rank[j]=copy;
	}
	if (rank==0)cout<<rank<<"\tv_rank mescolato:\n";
	for (int i=0; i<size; i++){
		if (rank==0)cout<<rank<<"\t"<<v_rank[i]<<endl;
	}
		
	MPI_Barrier(MPI_COMM_WORLD);

	for (int s=0; s<Ngen; s++){
		avL2=0;
		ordina(Pop, N, Ndim);
		if (s%Nmigr==0){
			for (int i=0; i<size-1; i+=2){
				for(int j=1; j<capitals.size(); j++){
					if (rank == v_rank[i]) {
						MPI_Send(&Pop[0][j][x], 1, MPI_DOUBLE, v_rank[i+1], x, MPI_COMM_WORLD);
						MPI_Send(&Pop[0][j][y], 1, MPI_DOUBLE, v_rank[i+1], y, MPI_COMM_WORLD);
					 }else if (rank == v_rank[i+1]){
						MPI_Recv(&Pop[Ndim-1][j][x], 1, MPI_DOUBLE, v_rank[i], x, MPI_COMM_WORLD, &status);
						MPI_Recv(&Pop[Ndim-1][j][y], 1, MPI_DOUBLE, v_rank[i], y, MPI_COMM_WORLD, &status);
					}
					MPI_Barrier(MPI_COMM_WORLD);
					
					if (rank == v_rank[i+1]) {
						MPI_Send(&Pop[0][j][x], 1, MPI_DOUBLE, v_rank[i], x, MPI_COMM_WORLD);
						MPI_Send(&Pop[0][j][y], 1, MPI_DOUBLE, v_rank[i], y, MPI_COMM_WORLD);
					 }else if (rank == v_rank[i]){
						MPI_Recv(&Pop[Ndim-1][j][x], 1, MPI_DOUBLE, v_rank[i+1], x, MPI_COMM_WORLD, &status);
						MPI_Recv(&Pop[Ndim-1][j][y], 1, MPI_DOUBLE, v_rank[i+1], y, MPI_COMM_WORLD, &status);
					}
					MPI_Barrier(MPI_COMM_WORLD);
				}
				CalculateL2(Pop[Ndim-1],N);
				MPI_Barrier(MPI_COMM_WORLD);
			}
			ordina(Pop, N, Ndim);
			MPI_Barrier(MPI_COMM_WORLD);
		}
		NewPop(Pop, N, Ndim, p_c, p_m, rnd);
		ordina(Pop, N, Ndim);
		double bestL2send=Pop[0][0][0];
		for(int j=0; j<Ndim/2; j++){
	  		avL2+=Pop[j][0][0]/((Ndim/2)*size);
	  	}
		double bestL2rec[10];
		double avL2rec[10];
		for (int i=0; i<10; i++){
			bestL2rec[i]=1000000.0;
			avL2rec[i]=1000000.0;
		}
		bestL2rec[0]=bestL2send;
		MPI_Gather(&bestL2send,1,MPI_DOUBLE,bestL2rec,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(&avL2,1,MPI_DOUBLE,avL2rec,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (rank==0){
			for(int j=1; j<size; j++){
				avL2+=avL2rec[j];
			}
			bestL2send=massimo(bestL2rec, 10);
			bestaveregeL2.open("averege_best_L2_"+to_string(size)+".dat",ios::app);
	    		bestaveregeL2 << s+1 << "\t" << bestL2send << endl;
	  		bestaveregeL2.close();
	  		averegeL2.open("averege_L2_"+to_string(size)+".dat",ios::app);
	  		averegeL2 << s+1 << "\t" << avL2 << endl;
	  		averegeL2.close();
			cout<<"GENERATION  "<<s<<endl;
		}
	}
	namefile="final_Path_";
	namefile=namefile+to_string(size);
	WriteBestPath(Pop, N, citypathfinal, namefile ,rank);

	MPI_Finalize();
	return 0;
}











































//HO FATTO UN CODICE CHE PRENDE LE CITTÀ IN ENTRATA E CI FA UNA LISTA UGUALE PER TUTTI I RANK E MESCOLA LORDINE DEL VETTORE RANK UGUALE PER TUTTI

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>
#include <cstdlib>
#include "lib.h"
#include "mpi.h"
using namespace std;

int x=0;
int y=1;
int L2=0;

	
int main(int argc, char* argv[]){

	int size, rank;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	int N=0;
	int Ngen=500;
	int Ndim=500;
	int Nmigr=20;
	double p_c=75./100.;
	double p_m=10/100.;
	double avL2=0;
	string namefile;

	ofstream citypathfinal, citypathstart, averegeL2, bestaveregeL2;

	Random rnd;				//random element unico per ogni processo (evoluzione dei sistemi separata e indipendente)
	Random rnd_com;				//random element comune a tutti i processi, serve a generare liste per scambio tra processi 
	
//Random rnd/rnd_com initialization
	int seed0[4];
	int p01, p02;
	ifstream Primes0("Primes");
	if (Primes0.is_open()){
		Primes0 >> p01 >> p02;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes0.close();
	
	ifstream input0("seed.in");
	string property0;
	if (input0.is_open()){
		while ( !input0.eof() ){
			input0 >> property0;
			if( property0 == "RANDOMSEED" ){
		         	input0 >> seed0[0] >> seed0[1] >> seed0[2] >> seed0[3];
		         	rnd_com.SetRandom(seed0,p01,p02);
		      	}
		}
		input0.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	for(int i=0; i<size; i++){
		if(rank == i){
			ifstream input("seed_" + to_string(i) + ".in");
			string property;
			if (input.is_open()){
	   			while ( !input.eof() ){
	      			input >> property;
	      			if( property == "RANDOMSEED" ){
	        			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	        			rnd.SetRandom(seed,p1,p2);
	      			}
	   			}
	  			input.close();
			} else cerr << "PROBLEM: Unable to open seed.in" << endl;
			
		}
	}
//fine inizializzazione random

//lettura da file
	ifstream ifile;
	ifile.open("American_capitals.in");
	vector<double> city{0.,0.};
	vector<double>initialL2{0.};
	vector<vector<double>> capitals;
	capitals.push_back(initialL2);
	//check to see that the file was opened correctly:
	if (!ifile.is_open()) {
		cerr << "There was a problem opening the input file: American_capitals.in!\n";
		exit(1);
	}
	double longitude = 0.0;
	double latitude = 0.0;
	string word;
	
	ifile >> word>> word>> word>> word;
	while (ifile >> word >> word >> longitude >> latitude) {
		city[x]=longitude;
		city[y]=latitude;
		capitals.push_back(city);
	}

	N=capitals.size()-1;

	MPI_Barrier(MPI_COMM_WORLD);
	vector<vector<vector<double>>> Pop;
	Pop = GeneraPopStart(capitals, N, Ndim, rnd);
	ordina(Pop, N, Ndim);
	namefile="start_Path_";
	namefile=namefile+to_string(size);
	WriteBestPath(Pop, N, citypathstart, namefile ,rank);
	vector<int> v_rank;
	
	if (rank==0)cout<<rank<<"\tv_rank:\n";
	for (int i=0; i<size; i++){
		v_rank.push_back(i);
		if (rank==0)cout<<rank<<"\t"<<v_rank[i]<<endl;
	}
	for (int i=0; i<size; i++){
		int j=(int)(rnd_com.Rannyu(0,size));
		int copy=v_rank[i];
		v_rank[i]=v_rank[j];
		v_rank[j]=copy;
	}
	if (rank==0)cout<<rank<<"\tv_rank mescolato:\n";
	for (int i=0; i<size; i++){
		if (rank==0)cout<<rank<<"\t"<<v_rank[i]<<endl;
	}
		
	MPI_Barrier(MPI_COMM_WORLD);
	double x_stock=0.;
	double y_stock=0.;
	for (int s=0; s<Ngen; s++){
		if (rank==0)cout<<s<<endl;
		avL2=0;
		ordina(Pop, N, Ndim);
		if (s%Nmigr==0&&s!=0&&Nmigr!=1){
			if (rank==0)cout<<"entro nel ciclo di scambio"<<endl;
			for (int i=0; i<size-1; i+=2){
				if (rank==0)cout<<"\t"<<i<<endl;
				for(int j=1; j<capitals.size(); j++){
					if (rank==0)cout<<"\t"<<"\t"<<j<<endl;
					x_stock=0.;
					y_stock=0.;
					if (rank == v_rank[i]) {
						MPI_Send(&Pop[0][j][x], 1, MPI_DOUBLE, v_rank[i+1], x, MPI_COMM_WORLD);
						MPI_Send(&Pop[0][j][y], 1, MPI_DOUBLE, v_rank[i+1], y, MPI_COMM_WORLD);
						
					 }else if (rank == v_rank[i+1]){
						MPI_Recv(&x_stock, 1, MPI_DOUBLE, v_rank[i], x, MPI_COMM_WORLD, &status);
						MPI_Recv(&y_stock, 1, MPI_DOUBLE, v_rank[i], y, MPI_COMM_WORLD, &status);
					}
					MPI_Barrier(MPI_COMM_WORLD);
					
					if (rank == v_rank[i+1]) {
						MPI_Send(&Pop[0][j][x], 1, MPI_DOUBLE, v_rank[i], x, MPI_COMM_WORLD);
						MPI_Send(&Pop[0][j][y], 1, MPI_DOUBLE, v_rank[i], y, MPI_COMM_WORLD);
					 }else if (rank == v_rank[i]){
						MPI_Recv(&x_stock, 1, MPI_DOUBLE, v_rank[i+1], x, MPI_COMM_WORLD, &status);
						MPI_Recv(&y_stock, 1, MPI_DOUBLE, v_rank[i+1], y, MPI_COMM_WORLD, &status);
					}
					Pop[0][j][x]=x_stock;
					Pop[0][j][y]=y_stock;
					MPI_Barrier(MPI_COMM_WORLD);
					if (rank==0)cout<<"\t"<<"\t"<<"\t"<<"sostituiti"<<endl;
				}
				CalculateL2(Pop[0],N);
				MPI_Barrier(MPI_COMM_WORLD);
			}
			ordina(Pop, N, Ndim);
			MPI_Barrier(MPI_COMM_WORLD);
		}
		NewPop(Pop, N, Ndim, p_c, p_m, rnd);
		ordina(Pop, N, Ndim);
		double bestL2send=Pop[0][0][0]/size;
		for(int j=0; j<Ndim; j++){
	  		avL2+=Pop[j][0][0]/((Ndim)*size);
	  	}
		double bestL2rec[10];
		double avL2rec[10];
		for (int i=0; i<10; i++){
			bestL2rec[i]=1000000.0;
			avL2rec[i]=1000000.0;
		}
		bestL2rec[0]=bestL2send;
		MPI_Gather(&bestL2send,1,MPI_DOUBLE,bestL2rec,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(&avL2,1,MPI_DOUBLE,avL2rec,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (rank==0){
			for(int j=1; j<size; j++){
				bestL2send+=bestL2rec[j];
				avL2+=avL2rec[j];
			}
			//bestL2send=massimo(bestL2rec, 10);
			bestaveregeL2.open("averege_best_L2_"+to_string(size)+".dat",ios::app);
	    		bestaveregeL2 << s+1 << "\t" << bestL2send << endl;
	  		bestaveregeL2.close();
	  		averegeL2.open("averege_L2_"+to_string(size)+".dat",ios::app);
	  		averegeL2 << s+1 << "\t" << avL2 << endl;
	  		averegeL2.close();
			cout<<"GENERATION  "<<s<<endl;
		}
	}
	ordina(Pop, N, Ndim);
	namefile="final_Path_";
	namefile=namefile+to_string(size);
	WriteBestPath(Pop, N, citypathfinal, namefile ,rank);

	MPI_Finalize();
	return 0;
}
