/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include <iomanip>


int main(){
    
    Input();             //Inizialization
 
    int nconf = 1;
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep) //questi sono gli step per blocco
        {
            Move();           //Move particles with Verlet algorithm
            if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
            if(istep%10 == 0){
                Measure();
                Accumulate(); //Update block averages
            //  ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                nconf += 1;
            }
        }
        Averages(iblk);   //Print results for current block
    }
    ConfFinal(); //Write final configuration
    
    return 0;
}

// ---------------------------------------------------------- //

double error (double AV, double AV2, int n){
    if (n==0) return 0;
    else return pow(((AV2 - AV*AV)/n), 0.5);
}

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf,ReadOldConf;
  //double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
    //ReadInput.open("input_solid.dat"); //Read input solid
    //ReadInput.open("input_liquid.dat"); //Read input liquid
    ReadInput.open("input_gas.dat"); //Read input gas
  
  ReadInput >> fase;
    if (fase==1) k = "solid";
    if (fase==2) k = "liquid";
    if (fase==3) k = "gas";
    
  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblk;
  ReadInput >> iprint;
  ReadInput >> irestart;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

   
//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

    //measurement of g(r)
    igofr = 4;
    nbins = 100;
    n_props = n_props + nbins;
    bin_size = (box/2.0)/(double)nbins;
    
    //Read initial configuration
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();
    
    if (irestart==0) {
           // fa l'algoritmo normale
        
        //Prepare initial velocities
        cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
        double sumv[3] = {0.0, 0.0, 0.0};
        for (int i=0; i<npart; ++i){
          vx[i] = rand()/double(RAND_MAX) - 0.5;
          vy[i] = rand()/double(RAND_MAX) - 0.5;
          vz[i] = rand()/double(RAND_MAX) - 0.5;

          sumv[0] += vx[i];
          sumv[1] += vy[i];
          sumv[2] += vz[i];
        }
        for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
        double sumv2 = 0.0, fs;
        for (int i=0; i<npart; ++i){
          vx[i] = vx[i] - sumv[0];
          vy[i] = vy[i] - sumv[1];
          vz[i] = vz[i] - sumv[2];

          sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        sumv2 /= (double)npart;

        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
        for (int i=0; i<npart; ++i){
          vx[i] *= fs;
          vy[i] *= fs;
          vz[i] *= fs;

          xold[i] = Pbc(x[i] - vx[i] * delta);
          yold[i] = Pbc(y[i] - vy[i] * delta);
          zold[i] = Pbc(z[i] - vz[i] * delta);
        }
       }
       else if (irestart==1){
           // prepara con le condizioni da old.0
           
           cout << "Read initial conditions from old.0" << endl << endl;
           ReadOldConf.open("old.0");
           for (int i=0; i<npart; ++i){
             ReadOldConf >> xold[i] >> yold[i] >> zold[i];
             xold[i] = xold[i] * box;
             yold[i] = yold[i] * box;
             zold[i] = zold[i] * box;
           }
           ReadOldConf.close();
           
           double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];
           for(int i=0; i<npart; ++i){ //Force acting on particle i
             fx[i] = Force(i,0);
             fy[i] = Force(i,1);
             fz[i] = Force(i,2);
           }

           double Temp, fs, sumv2 = 0.0;
           for(int i=0; i<npart; ++i){ //Verlet integration scheme

             xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
             ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
             znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

            vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
            vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
            vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

               x[i] = xnew;
               y[i] = ynew;
               z[i] = znew;
               
             sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
           }
           sumv2 /= (double)npart;
           
           Temp = sumv2/3.;
           fs = pow(temp/Temp, 0.5);
           
            for (int i=0; i<npart; i++){
                
               vx[i] *= fs;
               vy[i] *= fs;
               vz[i] *= fs;

               xold[i] = Pbc(x[i] - vx[i] * delta);
               yold[i] = Pbc(y[i] - vy[i] * delta);
               zold[i] = Pbc(z[i] - vz[i] * delta);
            }
       }


   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

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
  }
  return;
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

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
    
  ofstream Epot, Ekin, Etot, Temp;
    Epot.open(k+"/output_epot.dat",ios::app);
    Ekin.open(k+"/output_ekin.dat",ios::app);
    Temp.open(k+"/output_temp.dat",ios::app);
    Etot.open(k+"/output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
    
    //reset the hystogram of g(r)
    for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;
    
//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

        //update of the histogram of g(r)
        if (dr<box/2.)
        {
            bin = (int)(dr/bin_size); // poichè il bin che mi dà la distanza dr è quello tale che bin*bin_size = dr
            walker[bin + igofr] += 2;
        }
        
     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }        
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    walker[iv] = v;
    stima_pot = v/(double)npart; //Potential energy per particle
    Epot << stima_pot << endl;
    
    walker[ik] = t;
    stima_kin = t/(double)npart; //Kinetic energy per particle
    Ekin << stima_kin  << endl;
    
    walker[it]= (2.0 / 3.0) * t;
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    Temp << stima_temp << endl;

    walker[ie] = t+v;
    stima_tot = (t+v)/(double)npart; //Total energy per particle
    Etot << stima_tot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

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
    
   //double r, gdir;
   ofstream Gofr, Gave, Epot, Ekin, Etot, Temp;

   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    
    Epot.open(k+"/ave_epot.dat",ios::app);
    Ekin.open(k+"/ave_ekin.dat",ios::app);
    Etot.open(k+"/ave_etot.dat",ios::app);
    Temp.open(k+"/ave_temp.dat",ios::app);
    Gofr.open(k+"/ave_gofr.dat",ios::app);
    Gave.open(k+"/ave_gave.dat",ios::app);
    
    //Potential energy per particle
    stima_pot = blk_av[iv]/blk_norm/(double)npart; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;

    // Kinetic energy per particle
    stima_kin = blk_av[ik]/blk_norm/(double)npart; //Potential energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
    
    // Total energy per particle
    stima_tot = blk_av[ie]/blk_norm/(double)npart; //Potential energy
    glob_av[ie] += stima_tot;
    glob_av2[ie] += stima_tot*stima_tot;
    err_tot=Error(glob_av[ie],glob_av2[ie],iblk);
    Etot << setw(wd) << iblk <<  setw(wd) << stima_tot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_tot << endl;
    
    // Kinetic energy per particle
    stima_temp = blk_av[it]/blk_norm/(double)npart; //Potential energy
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;

//g(r)
    for (int i=0; i<nbins; i++){
        // r+dr è bin_size*(i+1) e r è i*bin_size
        // questo perchè per come ho costruito ho che (l/2)=nbins*bin_size che sarebbe il lato estremo dell'istogramma
        double deltaV = (4.*M_PI/3.)*(pow(bin_size*(i+1), 3.)-pow(bin_size*i, 3.));
        stima_gdir = blk_av[igofr+i]/blk_norm/(double)npart/rho/deltaV;
        glob_av[igofr + i] += stima_gdir;
        glob_av2[igofr + i] += stima_gdir*stima_gdir;
    }
    
    for (int i=0; i<nbins; i++){
        Gofr << setw(wd) << iblk << setw(wd) << (i*bin_size+bin_size/2.) << setw(wd) << glob_av[igofr + i]/(double)iblk << endl;
    }
    
    
    if(iblk==nblk) //Gave stampo solo l'ultimo blocco
    {
        for (int i=0; i<nbins; i++){
            err_gdir = Error(glob_av[igofr + i], glob_av2[igofr + i], iblk);
            Gave << setw(wd) << iblk << setw(wd) << (i*bin_size+bin_size/2.) << setw(wd) << glob_av[igofr + i]/(double)iblk << setw(wd) << err_gdir << endl;
        }
    }

    cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Gofr.close();
    Gave.close();
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
    
    cout << "Print preview configuration to file old.final " << endl << endl;
    WriteConf.open("old.final");

    for (int i=0; i<npart; ++i){
      WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
    }
    WriteConf.close();
  return;
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

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
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
