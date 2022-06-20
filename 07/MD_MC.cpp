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
#include "MD_MC.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  int nconf = 1;
//*
  for(int istep=1; istep <= 6000; istep++) 	//MD:serve per far partire la simulzione senza considerare il tempo di equilibrazione prima,2000 sol 6000 per liq e 120000 per gas
    {						//MC:serve per far partire la simulzione senza considerare il tempo di equilibrazione prima,1000 sol 2000 per liq e 30000 per gas
      Move();
    }
//*/
  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
      if(istep%10 == 0){
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  //cout << "Classic Lennard-Jones fluid        " << endl;
  //cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
  //cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  //cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  //cout << "The program uses Lennard-Jones units " << endl;

//Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

//Read input informations
  ReadInput.open("input.in");

  ReadInput >> iNVET;
  ReadInput >> restart;

  if(restart) Seed.open("seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  ReadInput >> temp;
  beta = 1.0/temp;
  //cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  //cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  //cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);  
  //cout << "Volume of the simulation box = " << vol << endl;
  //cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  //cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
  vtail=(8*pi*rho)*(1/(9*pow(rcut,9))-1/(3*pow(rcut,3)));
  ptail=(32*pi*rho)*(1/(9*pow(rcut,9))-1/(6*pow(rcut,3)));//calcolo delle tail correction
  ptail=ptail/vol;
    
  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  //cout << "The program perform Metropolis moves with uniform translations" << endl;
  //cout << "Moves parameter = " << delta << endl;
  //cout << "Number of blocks = " << nblk << endl;
  //cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements
  iv = 0; //Potential energy
  it = 1; //Temperature
  ik = 2; //Kinetic energy
  ie = 3; //Total energy
  iw = 4;//Pressure////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  n_props = 5; //Number of observables/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//measurement of g(r)
  ig = 5; //ig serve per il walker [] per partire da dopo le grandezze iv it ik ie iw
//  n_props = n_props + 1; //aggiorno n_props perché ho aggiunto il ig
  nbins = 100;
  n_props = n_props + nbins;//riaggiorno n_props perché ig in realtà è un mega vettorone da ig in poi
  bin_size = (box/2.0)/(double)nbins;
  
//Read initial configuration
  cout << "Read initial configuration" << endl << endl;
  if(restart)
  {
    ReadConf.open("config.out");
    ReadVelocity.open("velocity.out");
    for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
  }
  else 
  {
    ReadConf.open("config.in");
    cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i)
    {
      vx[i] = rnd.Gauss(0.,sqrt(temp));
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i)
    {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    cout << "velocity scale factor: " << fs << endl << endl;
    for (int i=0; i<npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();

  for (int i=0; i<npart; ++i)
  {
    if(iNVET)
    {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else
    {
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  
//Evaluate properties of the initial configuration
  Measure();

//Print initial values for measured properties
  //cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  //cout << "Initial temperature      = " << walker[it] << endl;
  //cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << endl;
  //cout << "Initial total energy     = " << walker[ie]/(double)npart << endl;
  //cout << "Initial pressure	    = " << walker[iw] << endl;/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  return;
}


void Move()
{
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) // Monte Carlo (NVT) move
  {
    for(int i=0; i<npart; ++i)
    {
    //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
      o = (int)(rnd.Rannyu()*npart);

    //Old
      energy_old = Boltzmann(x[o],y[o],z[o],o);

    //New
      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

      energy_new = Boltzmann(x[o],y[o],z[o],o);

    //Metropolis test
      p = exp(beta*(energy_old-energy_new));
      if(p >= rnd.Rannyu())  
      {
      //Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
        accepted = accepted + 1.0;
      } else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted = attempted + 1.0;
    }
  } else // Molecular Dynamics (NVE) move
  {
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i){ //Force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    for(int i=0; i<npart; ++i){ //Verlet integration scheme

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0;
      attempted = attempted + 1.0;
    }
  }
  return;
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
 
  return f;
}

void Measure() //Properties measurement
{
  int bin;
  double v = 0.0, kin=0.0, w=0.0;	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////, w=0.0;
  double vij, wij;	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////, wij;
  double dx, dy, dz, dr;
  
//reset the hystogram of g(r)
  for (int k=ig; k<ig+nbins; ++k) walker[k]=0.0;


//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {
// distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      
//update of the histogram of g(r)
      if( dr < box/2.0){
        bin = dr/bin_size;
        walker[bin+ig] = walker[bin+ig] + 2. ;
      }

      if(dr < rcut)
      {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        v += vij;
        wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);	
        w += wij;	
      }
    }          
  }

  for (int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  walker[iv] = 4.0 * v; // Potential energy
  walker[ik] = kin; // Kinetic energy
  walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature
  walker[ie] = 4.0 * v + kin;  // Total energy;
  walker[iw] = rho*walker[it]+(48.0/(3.0*vol))*w;	

  return;
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
   double gdir, norm, r;
   ofstream Gofr, Gave, Epot, Ekin, Etot, Temp, Pres;
   //const int wd=12;
    //epot_insta_GAS, epot_insta_LIQ, epot_insta_SOL
    Epot.open("output_epot_LIQ_MD.0",ios::app);
    //Ekin.open("output_ekin.dat",ios::app);
    //Temp.open("output_temp.dat",ios::app);
    //Etot.open("output_etot.dat",ios::app);
    Pres.open("output_pres_LIQ_MD.0",ios::app);
    Gofr.open("output_gofr_LIQ_MD.0",ios::app) ;
    Gave.open("output_gave_LIQ_MD.0",ios::app) ;
    
    if (iblk!=0){
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_pres = blk_av[iw]/blk_norm + ptail*(double)npart; ///////////////////////////////
    glob_av[iw] += stima_pres;    		   
    glob_av2[iw] += stima_pres*stima_pres;	   
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);
    
/*  stima_kin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = blk_av[ie]/blk_norm/(double)npart; //Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

    stima_temp = blk_av[it]/blk_norm; //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk); 
    
    
    //*/
//Potential energy per particle
    Epot << "\t" << iblk <<  "\t" << stima_pot << "\t" << glob_av[iv]/(double)iblk << "\t" << err_pot << endl;
//Pressure
    Pres << "\t" << iblk <<  "\t" << stima_pres << "\t" << glob_av[iw]/(double)iblk << "\t" << err_press << endl;
//Kinetic energy
/*    Ekin << "\t" << iblk <<  "\t" << stima_kin << "\t" << glob_av[ik]/(double)iblk << "\t" << err_kin << endl;
//Total energy
    Etot << "\t" << std::scientific << iblk <<  "\t" << stima_etot << "\t" << glob_av[ie]/(double)iblk << "\t"<< err_etot << endl;	////////////////////////////////////////
//Temperature
    Temp << "\t" << iblk <<  "\t" << stima_temp << "\t" << glob_av[it]/(double)iblk << "\t" << err_temp << endl;

    cout << "----------------------------" << endl << endl;
//*/


//g(r)
    for(int k=ig; k<n_props; k++){
      r = (k-ig)*bin_size;
      gdir = blk_av[k]/blk_norm;		
      glob_av[k] += gdir;
      glob_av2[k] += gdir*gdir;
      norm = 1./((4./3.)*pi*(pow((r + bin_size),3) - pow(r,3))*rho*npart);
      Gofr << norm*gdir << "  " ;
    }
    Gofr << endl;
		
		
    if(iblk == nblk){
      for(int k=ig; k<n_props; k++){
        r = (k-ig)*bin_size;
    	norm = 1./((4./3.)*pi*(pow((r + bin_size),3) - pow(r,3))*rho*npart);
    	Gave <<  r + bin_size/2.<< "\t" << norm*glob_av[k]/(double)nblk << "\t" << norm*Error(glob_av[k],glob_av2[k],nblk) << endl;
      }
    }
    }
    
    Epot.close();
    //Ekin.close();
    //Etot.close();
    //Temp.close();
    Pres.close();
    Gofr.close();
    Gave.close();
}


void ConfFinal(void)
{
  ofstream WriteConf, WriteVelocity, WriteSeed;

  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open("config.out");
  WriteVelocity.open("velocity.out");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}



/*  
  int main()
{ 
  Input(); //Inizialization
  int nconf = 1;
  //Averages(0);
  
  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    int istep=1;
    int eqstep=5000;
    for(istep; istep <= eqstep; istep++) 	//serve per far partire la simulzione senza considerare il tempo di equilibrazione prima
    {
      Move();
    }
    for(istep; istep <= nstep; istep++)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
      if(istep%10 == 0){
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}*/
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
