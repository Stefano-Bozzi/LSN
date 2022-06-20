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
#include "Monte_Carlo_QM_1D.h"

using namespace std;

int main(){
	cout<<endl;
  Input(); //Inizialization
  x_start=x_old;
  if (SA==1){
  	simulated_annealing();
  	ReadSAResult();
  }
  cout<<"starting simulation with best parameters..."<<endl;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(1,mu,sigma);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  rnd.SaveSeed();
  return 0;
}

  
  

/////FUNCTIONS//////

void Input(void)
{
//cout <<"starting input"<<endl;
  ifstream ReadInput, Primes, Seed;

//Read seed for random numbers
   int p1, p2;
   Primes.open("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   Seed.open("seed.in");
   Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   Seed.close();
  
//Read input informations
  ReadInput.open("input.dat");

	ReadInput >> SA;
  cout<<"simulation with simulated annealing? SA=1=>yes, SA=0=>no::: SA= "<<SA<<endl;

  ReadInput >> x_old;
  cout<<"x_start: "<<x_old<<endl;
  
  ReadInput >> sigma;
  cout<<"sigma start: "<<sigma<<endl;
  
  ReadInput >> mu;
  cout<<"mu start: "<<mu<<endl;
  
  ReadInput >> h;
  cout<<"h_tagliato: "<<h<<endl;
  
  ReadInput >> m;
  cout<<"m: "<<m<<endl;  

  ReadInput >> delta;
  cout<<"delta x: "<<delta<<endl;
  
  ReadInput >> nblk;
  cout<<"nblk: "<<nblk<<endl;

  ReadInput >> nstep;
  cout<<"nstep: "<<nstep<<endl;
  
/* 
  ReadInput >> beta_start;
  cout<<"beta_start: "<<beta_start<<endl;
  
  ReadInput >> beta_fin;
  cout<<"beta_fin: "<<beta_fin<<endl;
  
//*/  

	ReadInput >> temp;
  cout<<"temp: "<<temp<<endl;
  
  ReadInput >> delta_temp;
  cout<<"delta_temp: "<<delta_temp<<endl;
  
  ReadInput >> ntemp;
  cout<<"num_sim_per_temp: "<<ntemp<<endl;
  
  ReadInput >> nstep_SA;
  cout<<"nstep_SA: "<<nstep_SA<<endl;
  
  ReadInput >> delta_mu;
  cout<<"delta_mu: "<<delta_mu<<endl;
  
  ReadInput >> delta_sigma;
  cout<<"delta_sigma: "<<delta_sigma<<endl;
  
  ReadInput.close();


//*
//Prepare arrays for measurements
  ih = 0; //Hamiltonian
  
  n_props = 1; //Number of observables
  
//Evaluate Hamiltonian etc. of the initial configuration
  Measure();
 
	cout<<"walker_ih start from: "<<walker[ih]<<endl<<endl;
}

double Psi(double x, double mu, double sigma){
//exp(-pow( (x-x0)/sigma, 2.) /2.)
  return exp(-(pow(x-mu,2))/(2*pow(sigma,2))) + exp(-(pow(x+mu,2))/(2*pow(sigma,2)));
}

double Prob(double x, double mu, double sigma)
{
  double P =pow(Psi(x,mu,sigma),2);
  return P;
}

double secondDeriv(double x, double mu, double sigma){
	return -(1/pow(sigma,2))*(Psi(x,mu,sigma)-(1/pow(sigma,2))*(pow(x-mu,2)*exp(-(pow(x-mu,2))/(2*pow(sigma,2))) + pow(x+mu,2)*exp(-(pow(x+mu,2))/(2*pow(sigma,2)))));
}
double Pot(double x){
  double v = pow(x,4)-(5./2.)*pow(x,2);
  //double v = 0.5*x*x;
  return v;
}

void Move(int a, double mu, double sigma){
  double x_new;
  attempted++;
  x_new=rnd.Rannyu(x_old-delta,x_old+delta);
  double Prob_old=Prob(x_old, mu, sigma);
  double Prob_new=Prob(x_new, mu, sigma);
  double alfa=min(1.,Prob_new/Prob_old);
  if( rnd.Rannyu() < alfa ){
    x_old=x_new;
    accepted++;
  }
  if (a==1){
    ofstream H;
    if (SA==0) H.open("positions.out",ios::app);
    else if (SA==1) H.open("positions_SA.out",ios::app);
    H<<std::setprecision(9)<<x_old<<endl;
    H.close(); 
  }
}

void Measure(){
  walker[ih]=((-pow(h,2)/(2*m))*(secondDeriv(x_old,mu,sigma)/Psi(x_old,mu,sigma))+Pot(x_old));
}

void simulated_annealing(void){
	cout<<"starting simulated annealing ..."<<endl;
	double H_old=0.0;
	double H_new=0.0;
	double H_old_2=0.0;
	double H_new_2=0.0;
	double c=temp/delta_temp;
	double j=0;
	x_old=x_start;
	//for(double beta=beta_start; beta<=beta_fin; beta+=0.5){
	for(double t=temp; t>=delta_temp; t-=delta_temp){
		//temp=1/beta;
		//cout<<std::setprecision(2)<<(100*t)/delta_temp<<"%"<<endl;
		double beta=1/t;
		for (int i=0; i<ntemp; i++){
			H_old=0.0;
			H_new=0.0;
			H_old_2=0.0;
			H_new_2=0.0;
			x_old=x_start;
		  for (int n=0; n<nstep_SA; n++){
		    Move(0,mu,sigma);
		    H_old+=((-pow(h,2)/(2*m))*(secondDeriv(x_old,mu,sigma)/Psi(x_old,mu,sigma))+Pot(x_old));
		    H_old_2+=((-pow(h,2)/(2*m))*(secondDeriv(x_old,mu,sigma)/Psi(x_old,mu,sigma))+Pot(x_old))*((-pow(h,2)/(2*m))*(secondDeriv(x_old,mu,sigma)/Psi(x_old,mu,sigma))+Pot(x_old));
		  }
		  H_old=H_old/nstep_SA;
		  H_old_2=H_old_2/nstep_SA;
				mu_new =fabs( mu + delta_mu*(rnd.Rannyu() -0.5) );
        sigma_new =fabs( sigma + delta_sigma*(rnd.Rannyu() - 0.5) );
				
		  x_old=x_start;
		  for (int n=0; n<nstep_SA; n++){
		    Move(0,mu_new,sigma_new);
		    H_new+=((-pow(h,2)/(2*m))*(secondDeriv(x_old,mu_new,sigma_new)/Psi(x_old,mu_new,sigma_new))+Pot(x_old));
		    H_new_2+=((-pow(h,2)/(2*m))*(secondDeriv(x_old,mu_new,sigma_new)/Psi(x_old,mu_new,sigma_new))+Pot(x_old))*((-pow(h,2)/(2*m))*(secondDeriv(x_old,mu_new,sigma_new)/Psi(x_old,mu_new,sigma_new))+Pot(x_old));
		  }
		  H_new=H_new/nstep_SA;
		  H_new_2=H_new_2/nstep_SA;
		  double p=1.0;
		  if (H_new > H_old) p=exp(-beta*(H_new-H_old));
		  if(p >= rnd.Rannyu() ){
		  	mu = mu_new;
		  	sigma = sigma_new;
		  	H_old=H_new;
		  	H_old_2=H_new_2;
		    SA_accepted++;
		  }
		  SA_attempted++;
		}
		cout<<std::fixed;
		cout<<std::setprecision(1)<<100*j/c<<"%\tactual temp= "<<std::setprecision(4)<<t<<endl;
		j++;
		ofstream fout;
		fout.open("H_temp.out",ios::app);
		fout<<std::setprecision(9)<<t<<"\t"<<mu<<"\t"<<sigma<<"\t"<<H_old<<"\t"<<Error(H_old*nstep_SA, H_old_2*nstep_SA, nstep_SA)<<endl;
		fout.close();
	}
	cout<<"end of simulated annealing"<<endl;
}

void ReadSAResult(void){
	double t,H,H_best,H_err_best,H_err,mu_best,sigma_best;
	ifstream ReadIn;
	ReadIn.open("H_temp.out");
	if (!ReadIn.is_open()) {
		cerr << "There was a problem opening the input file: H_temp.dat!\n";
		exit(1);
	}
	ReadIn >> t >> mu_best >> sigma_best >> H_best >> H_err_best;
	cout<<"sorting parameters...";
	while (ReadIn >> t >> mu >> sigma >> H >> H_err) {
		if(H<H_best){
			H_best=H;
			mu_best=mu;
			sigma_best=sigma;
		}
	}
	sigma=sigma_best;
	mu=mu_best;
	cout<<"best parameters found"<<endl;
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


void Averages(int iblk) //Print results for current block
{
    
    ofstream Ham;
   
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    if (SA==0) Ham.open("output_Ham.out",ios::app);
    else if (SA==1) Ham.open("output_Ham_SA.out",ios::app);
    stima_h = blk_av[ih]/blk_norm; //Energy
    glob_av[ih]  += stima_h;
    glob_av2[ih] += stima_h*stima_h;
    err_h=Error(glob_av[ih],glob_av2[ih],iblk);
    Ham <<std::setprecision(9)<<iblk <<  "\t" << stima_h << "\t" << glob_av[ih]/(double)iblk << "\t" << err_h << endl;
    Ham.close();

    cout << "----------------------------" << endl << endl;
}

double min(double a, double b){
    if (a < b){
        return a;
    } else {
        return b;
    }
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
