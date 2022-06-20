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

////////////////////////////////////////////////////////////////////////////////LEGENDA DEL CODICE:////////////////////////////////////////////////////////////////////////
/*
ogni città è un vettore con coordinate [x,y]: vector<double>

cityCircle/Square: sono liste di città ordinate, fatte nel seguente modo: nella posizione 0 (chiamata anche L2) vi è la quantità L2 appunto, caratteristica della lista

									vector <vector <double>>: posizione  |  elemento  |  sottoelemento
												      0=L2       0 L2          1 none
													
													1        città 		0 [x|
																1 |y]
																
													2        città 		0 [x|
																1 |y]
													...
Popolazione: insieme di liste di città che formano una popolazione, tipologia:
									
								vector< vector <vector <double>>> : posizione  |     elemento      |  sottoelemento  |  sottoelemento
													0        cityCircle/Square         L2
													
																          città 	     [x|
													    					             |y]

													
																      	  città 	     [x|
													    						     |y]
													    			           ...
													1        cityCircle/Square         L2
													
																          città 	     [x|
													    					             |y]

													
																      	  città 	     [x|
													    						     |y]
													    			           ...
													    			           
ESEMPIO DI POIZIONE DI ELEMENTI:	elemento L2 di una lista cityCircle: cityCircle[0][0] (oppure cityCircle[L2][0])
					
					coordinata x della prima città di una lista cityCircle: cityCircle[1][x]
					
					coordinata y della seconda città di una lista cityCircle: cityCircle[2][y]
					
					elemento L2 della seconda lista cityCircle di una popolazione Pop: Pop[1][L2][0]
					
					coordinata y della prima città della terza lista cityCircle di una popolazione Pop: Pop[3][1][y]
//*/
					
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//MAIN PROGRAM_______________________________________________________________________________________________________________________________________________________________________________________


int main (int argc, char *argv[]){
	Random rnd;				//creating a random element to initialize the random generator
	rand_initialize(rnd);	
	int N=34;
	int Ngen=600;
	int Ndim=500;
	double p_c=75./100.;
	double p_m=20/100.;
	double avL2=0;
	ofstream best_L2, av_L2, cityPath, cityPath1, cityPathfinal;
	
	vector <vector<double>> citySquare, cityCircle;//, citySquareStart, cityCircleStart;
	vector <double> pos,initialL2;
	pos={0.,0.};
	initialL2={0.};
	rnd.Rannyu();
	citySquare.push_back(initialL2);
	cityCircle.push_back(initialL2);
	
	for (int i=0; i<N; i++){
		double theta=rnd.Rannyu(0.,2*M_PI);
		pos[x]=rnd.Rannyu();
		pos[y]=rnd.Rannyu();
//		pos[x]=i+1;
//		pos[y]=i+1;
		citySquare.push_back(pos);
		pos[x]=cos(theta);
		pos[y]=sin(theta);
		cityCircle.push_back(pos);	
	}
	
	CheckCoincidence(cityCircle,N,rnd);
	CheckCoincidence(citySquare,N,rnd);
	CalculateL2(cityCircle,N);
	CalculateL2(citySquare,N);
	vector<vector<vector<double>>> Pop;	
	//Pop = GeneraPopStart(cityCircle, N ,Ndim, rnd);
	Pop = GeneraPopStart(citySquare, N ,Ndim, rnd);
	
	cout<<"Writing first casual path for Square..."<<endl<<endl;
	WriteBestPath(Pop, N, cityPath1, "Square_first_Path" ,0);
	
	ordina(Pop, N, Ndim);
	for (int i=0; i<Ngen; i++){
		avL2=0;
		ordina(Pop, N, Ndim);
		if ((i==0||(i+1)%5==0)&&i<520){
			WriteBestPath(Pop, N, cityPath, "data_Square/cityPath_" ,i+1);
		}
		NewPop(Pop, N, Ndim, p_c, p_m, rnd);
		ordina(Pop, N, Ndim);
		//cout<<"Writing BEST L2 for Circle..."<<endl<<endl;
   		best_L2.open("best_L2_Square.dat",ios::app);
    		best_L2 << i+1 << "\t" << Pop[0][L2][0] << endl;
  		best_L2.close();
  		//*
  		cout<<"Writing AVEREGE L2 for Square..."<<endl<<endl;
  		av_L2.open("av_L2_Square.dat",ios::app);
  		for(int j=0; j<Ndim/2; j++){
  			avL2+=Pop[j][0][0]/(Ndim/2);
  		}
    		av_L2 << i+1 << "\t" << avL2 << endl;
  		av_L2.close();
  		//*/
  		cout<<"GENERATION  "<<i<<endl;
	}
	
	WriteBestPath(Pop, N, cityPathfinal, "Square_final_Path" ,0);
	
	return 0;
}

