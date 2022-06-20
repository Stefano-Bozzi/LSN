/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __QM__
#define __QM__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props,ih;
//double nbins;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_h;
double err_h;

//configuration
double sigma, mu, h, m;

// simulation
int nstep, nblk;
double delta, x_old;

//simulated annealing
double temp,delta_temp,delta_mu,delta_sigma;
int ntemp,nstep_SA,SA;
double sigma_new, mu_new, x_start;
double SA_accepted,SA_attempted;

//functions
void Input(void);
double Psi(double, double, double);
double Prob(double, double, double);
double secondDeriv(double, double, double);
double Pot(double);
void Move(int, double, double);
void Measure(void);
void simulated_annealing(void);
void ReadSAResult(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
double min(double, double);
double Error(double,double,int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
