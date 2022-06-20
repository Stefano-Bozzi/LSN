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
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(int argc, char** argv){


 if (argc!=2){cout<<"proper usage of the program:"<<endl<<argv[0]<<" 0\tif you want single T simulation"<<endl<<argv[0]<<" 1\tif you want all T simulation"<<endl;}	//SETUP E CONTROLLO
 int numtemp=atoi(argv[1]);
 int equ_metro=200;
 int equ_gibbs=300;
 int equ=0;
 if (metro==1){equ=equ_metro;}
 else if (metro!=0){equ=equ_gibbs;}
 
 
 if (numtemp==0){					//SE VOGLIO FARE UNA SIMULAZIONE PER SOLO UNA TEMP. PER ESEMPIO TEMPO EQ
  Input(); //Inizialization
  for (int i=0; i<=equ; i++){	//tempo di equilibrazione
     Move(metro,numtemp);
  }
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro,numtemp);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk,numtemp);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
 }
 
 
 else if(numtemp==1){					//SE VOGLIO SIMULAZIONE COMPLETA PER 11 TEMPERATURE
   
   for(int i=0; i<11; i++){  //ciclo sulle temperature
    Input(); //Inizialization
    temp=0.5+i*0.15;
    beta=1/temp;
    
    if (temp!=0.5){ //se la temperatura non è 00.5 allora devo ripartire da configurazione precedente
     equ_metro=50;
     equ_gibbs=20;
     equ=0;
     if (metro==1){equ=equ_metro;}
     else if (metro!=0){equ=equ_gibbs;}
     ifstream ReadInput;
     ReadInput.open("config.final");
     for (int i=0; i<nspin; ++i){
      ReadInput >> s[i];
     }
     ReadInput.close();
    }
    
    for (int i=0; i<=equ; i++){	//tempo di equilibrazione
     Move(metro,numtemp);
    }
    
    //simulazione:
    
    for(int iblk=1; iblk <= nblk; ++iblk){
     Reset(iblk);   //Reset block averages
     for(int istep=1; istep <= nstep; ++istep){
      Move(metro,numtemp);
      Measure();
      Accumulate(); //Update block averages
     }
     Averages(iblk, numtemp);   //Print results for current block
    }
    ConfFinal(); //Write final configuration
    
    ofstream FINAL;	//creo un file incolonnato correttamente per lettura finale rapida
    if (metro==1){FINAL.open("output.FINAL_0_02_Metro",ios::app);}
    else if (metro!=1){FINAL.open("output.FINAL_0_02_Gibbs",ios::app);}
    cout<<"Temp simulation done: "<<temp<<endl;
    FINAL << "\t" << temp <<  "\t" << glob_av[iu]/(double)nblk << "\t" << err_u << "\t" << glob_av[ic]/(double)nblk << "\t" << err_c << "\t" << glob_av[im]/(double)nblk << "\t" << err_m << "\t" << glob_av[ix]/(double)nblk << "\t" << err_x << endl;
    FINAL.close();
    cout<<"Temp simulation done: "<<temp<<endl;
   }
 }
  
  return 0;
}


/////FUNCTIONS//////

void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> restart;		
  beta = 1.0/temp;
  cout << "Restart = " << restart << endl;
  
  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  if (restart==0){
   for (int i=0; i<nspin; ++i)
   {
     if(rnd.Rannyu() >= 0.5) s[i] = 1;
     else s[i] = -1;
   }
  }
  else if(restart==1){
   ReadInput.open("config.final");
   for (int i=0; i<nspin; ++i)
   {
    ReadInput >> s[i];
   }
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro, int numtemp)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1){ //Metropolis
     attempted++;
     energy_old=Boltzmann(s[o], o);
     energy_new=Boltzmann(-s[o],o);
     double d_E=energy_new - energy_old;
     if( d_E < 0 ){s[o]=-s[o];accepted++;}
     else if(d_E>=0&&rnd.Rannyu() < exp(-beta*d_E)){
      s[o]=-s[o];
      accepted++;
     }
//*
     if (numtemp==0){
      Measure();
      Accumulate();
      ofstream EneM;
      EneM.open("output.eneM.0",ios::app);
      stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
      EneM <<std::scientific<<i <<  "\t" << stima_u << endl;
      EneM.close();
     }
//*/
    }
    else{
     attempted++;
     accepted++;
     energy_up = Boltzmann(1, o);
     energy_down = Boltzmann(-1, o);
     if(rnd.Rannyu() < 1/(1+exp(beta*(energy_down - energy_up))) ){
      s[o] = -1;
     }
     else{
      s[o] = 1;
     }
//*
     if (numtemp==0){
      Measure();
      Accumulate();
      ofstream EneG;
      EneG.open("output.eneG.0",ios::app);
      stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
      EneG <<std::scientific<<i <<  "\t" << stima_u << endl;
      EneG.close();
     }
//*/
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;		
  walker[im] = m;
  walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk, int numtemp) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   
    
    if(numtemp==0){
     cout << "Block number " << iblk << endl;
     cout << "Acceptance rate " << accepted/attempted << endl << endl;
    }
//*    
    Ene.open("outputG.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene <<std::scientific<<iblk <<  "\t" << stima_u << "\t" << glob_av[iu]/(double)iblk << "\t" << err_u << endl;
    Ene.close();
//*
    Heat.open("outputG.heat.0", ios::app);
    stima_c = beta*beta*(blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm,2))/(double)nspin;
    glob_av[ic] += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat <<std::scientific<< iblk <<  "\t" << stima_c << "\t" << glob_av[ic]/(double)iblk << "\t" << err_c << endl;
    Heat.close();
    
    Mag.open("outputG.mag.0", ios::app);
    stima_m = (blk_av[im]/blk_norm)/(double)nspin;
    glob_av[im] += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag <<std::scientific<< iblk <<  "\t" << stima_m << "\t" << glob_av[im]/(double)iblk << "\t" << err_m << endl;
    Mag.close();
    
    Chi.open("outputG.chi.0", ios::app);
    if (h!= 0.){stima_x = beta*(blk_av[ix]/blk_norm - pow(blk_av[im]/blk_norm,2))/(double)nspin;}
    else if (h==0.){stima_x = beta*(blk_av[ix]/blk_norm)/(double)nspin;}
    glob_av[ix] += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi <<std::scientific<< iblk <<  "\t" << stima_x << "\t" << glob_av[ix]/(double)iblk << "\t" << err_x << endl;
    Chi.close();
//*/

    if(numtemp==0){						/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     cout << "----------------------------" << endl << endl;
    }
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
