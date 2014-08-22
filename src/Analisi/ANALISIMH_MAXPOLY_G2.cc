//==================================================================
// COPYRIGHT 2013 Claudio Bonati
// e-mail: claudio.bonati82@gmail.com
//==================================================================
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//====================================================================


#include<cmath>
#include<cstdlib>
#include<ctime> 
#include<fstream>
#include<iostream>
#include<new>
#include<string>

#include"./random.cc"

#define REAL long double


/* global variables  to be read from inputfile */
std::string *dati_in;
std::string dati_out;

long int *block;
long int *sample;
long int *numblocks;
long int numblockstot;
REAL *betas;

REAL betamin;
REAL betamax;
REAL betastep;
REAL volume;

int numfiles;
long int campione;

/* const global variables  */

const int numobs=2;
const int numobs_tot_jack=4;

const REAL pi=3.141592653589793238462643383279502884197169399375105820974944;

const REAL energy_off=0.0;

/*-----*/

REAL logsum(REAL a, REAL b);
REAL logdiff(REAL a, REAL b);

/*--------*/


int readinput(char *in_file)
    {
    std::ifstream input;
    std::string str, str2;
    long int temp_li;
    int temp_i;
    REAL temp_re;
    const int char_to_ignore=500;

    input.open(in_file, std::ios::in);
    if(input.is_open())
      {
      for(int i=0; i<6; i++)
         {
         input >> str;

         if(str=="outfile")
           { 
           input >> str2;
           dati_out=str2;
           }
         else if(str=="betamin")
                {
                input >> temp_re;
                betamin=temp_re;
                }
         else if(str=="betamax")
                {
                input >> temp_re;
                betamax=temp_re;
                }
         else if(str=="betastep") 
                {
                input >> temp_re;
                betastep=temp_re;
                }
         else if(str=="volume") 
                {
                input >> temp_li;
                volume=(REAL) temp_li;
                }
         else if(str=="numfiles") 
                {
                input >> temp_i;
                numfiles=temp_i;
                }
         else
             {
             std::cerr << "Error: unrecognized option \"" << str << "\" in the file \"" << in_file << "\"\n";
             return 1;     
             }

         input.ignore(char_to_ignore, '\n');
         }

      dati_in = new std::string[numfiles];
      block = new long int [numfiles];
      sample = new long int [numfiles];
      numblocks = new long int [numfiles];
      betas = new REAL [numfiles];

      campione=0;
      numblockstot=0;
      for(int i=0; i<numfiles; i++)
         {
         input >> str2;
         dati_in[i]=str2;

         input >> temp_re;
         betas[i]=temp_re;
     
         input >> temp_li;
         sample[i]=temp_li;

         input >> temp_li;
         block[i]=temp_li;
   
         numblocks[i]=sample[i]/block[i]; 
         numblockstot+=numblocks[i];
         campione+=block[i]*numblocks[i];

         input.ignore(char_to_ignore, '\n');
         }

      input.close();

      return 0;
      }
    else
      {
      std::cerr << "Error in opening the file \"" << in_file <<"\"\n";
      return 1;
      }
    }



class Data {
private:
  REAL*** dati0; // complete data
  REAL** dati;  // bootstrapped data
  long int * sample_loc;  // bootstrapped equivalent of "sample"
  REAL*  medie;
  REAL* zetas0;
  REAL* zetas;
  long int campione_loc;

public:
  Data();
  ~Data();
  void initfromfile(std::string *nome_file);
  void bootstrapdata(void);
  void computezetas0(void);
  void computezetas(void);
  void computeobs(REAL beta);
  REAL betapc(REAL beta_low, REAL beta_high, REAL precision);
}; 



Data::Data(void)
  {
  long int i, j;
  
  dati = new REAL * [campione];
  for(i=0; i<campione; i++)
     {
     dati[i]=new REAL [numobs+1];  // in the +1 the beta value is stored
     }

  dati0 = new REAL ** [numfiles];
  for(i=0; i<numfiles; i++)
     {
     dati0[i] = new REAL * [sample[i]];
     for(j=0; j<sample[i]; j++)
        {
        dati0[i][j]=new REAL [numobs];
        }
     }

  sample_loc = new long int [numfiles];
 
  medie  = new REAL [numobs_tot_jack];

  zetas = new REAL [numfiles];
  zetas0 = new REAL [numfiles];
  }


Data::~Data(void)
  {
  long int i, j;
  
  for(i=0; i<campione; i++)
     {
     delete [] dati[i];
     }
  delete [] dati;

  for(i=0; i<numfiles; i++)
     {
     for(j=0; j<sample[i]; j++)
        {
        delete [] dati0[i][j];
        }
     delete [] dati0[i];
     }
  delete [] dati0;

  delete [] sample_loc;

  delete [] medie;

  delete [] zetas;
  delete [] zetas0;
  }



void Data::initfromfile(std::string *nome_file)
  {
  std::ifstream filein;
  int j;
  REAL temp1;
  long int i, f;

  std::cout << "\nReading data ... \n";

  for(f=0; f<numfiles; f++)
     {
     filein.open( (nome_file[f]).c_str(), std::ios::in);
     if(filein.is_open())
       {
       for(i=0; i<sample[f]; i++)
          {
          for(j=0; j<numobs; j++)
             {
             filein >> temp1;
             if(j==0)
               {
               dati0[f][i][j]=temp1 + energy_off; // so that it is >0
               }
             else
               {
               dati0[f][i][j]=temp1;
               }
             }
          }
       filein.close();
       }
     else
       {
       std::cerr << "Error in opening the file "<<nome_file[f]<<" \n";
       exit(1);
       }
     }

  }


void Data::bootstrapdata(void)
  {
  long int f, i, r, j;
  int stop, test;
  
  campione_loc=0;
  for(f=0; f<numfiles; f++)
     {
     sample_loc[f]=0;
     }

  stop=0;
  while(stop==0)
       {
       f=(long int)(numfiles*casuale());

       if(campione_loc+block[f]<campione)
         {
         r=(long int) (numblocks[f]*casuale());
         for(i=0; i<block[f]; i++)    
            {                                    
            for(j=0; j<numobs; j++)
               {
               dati[campione_loc+i][j]=dati0[f][r*block[f]+i][j];
               }
            dati[campione_loc+i][numobs]=betas[f];
            }
        
         campione_loc+=block[f];
         sample_loc[f]+=block[f];
         }
       else
         {
         stop=1;
         }
       }

  test=1;
  for(f=0; f<numfiles; f++)
     {
     if(sample_loc[f]==0)
       {
       test=0;
       }
     }

  if(test==0)
    {
    bootstrapdata();
    }
  }



// compute the partition functions of the global data
void Data::computezetas0(void) 
  {
  const REAL minvalue=1.0e-7;
  const REAL minvaluespeedup=1.0e-4;

  REAL stopcond, stopcondold, energy, temp1, temp2;
  int i, k, f, j, l, first, fileexists;
  REAL *zetasold;
  int speedup;
  long int count=1;

  std::ifstream infile;
  std::ofstream outfile;

  std::cout << "\nSelfconsistent computation started\n";

  zetasold = new REAL [numfiles];

  infile.open("zetas.dat", std::ios::in);
  if(infile.is_open())
    {
    std::cout << "reading stored zetas" << std::endl;
    for(i=0; i<numfiles; i++)
        {
        infile >> zetasold[i];
        }
    infile.close();
    speedup=1;
    fileexists=0;
    }
  else
    {
    for(i=0; i<numfiles; i++)
       {
       zetasold[i]=1.0+casuale();
       }
    speedup=500;
    fileexists=1;
    }

  std::cout << "remainder (stopvalue="<<minvalue<<"):" << std::endl;

  stopcond=1.0;
  while(stopcond>minvalue)
       {
       for(k=0; k<numfiles; k++) // estimate zetas[k]
          {
          first=0;

          for(f=0; f<numfiles; f++)               //
             {                                    // 
             for(j=0; j<sample[f]; j+=speedup)    //  this means for all the points 
                {                                 //
                           
                if(j+speedup<sample[f])
                  {     
                  energy=0.0;
                  for(l=0; l<speedup; l++)
                     {
                     energy+=7.0*volume*6.0*(1.0-dati0[f][j+l][0]);
                     }
                  energy/=speedup;
                  }
                else
                  {
                  energy=7.0*volume*6.0*(1.0-dati0[f][j][0]);
                  }
                temp1=log(sample[0]/speedup)-zetasold[0]+(betas[k]-betas[0])*energy;
                for(l=1; l<numfiles; l++)
                   {
                   temp2=log(sample[l]/speedup)-zetasold[l]+(betas[k]-betas[l])*energy;
                   temp1=logsum(temp1, temp2);
                   }

                if(first==0)
                  {
                  zetas0[k]=-temp1;
                  } 
                else
                  {
                  zetas0[k]=logsum(zetas0[k], -temp1);
                  }
                
                first=1;
                }
             } 
          }
      
       stopcondold=stopcond;
       stopcond=0.0;
       for(i=0; i<numfiles; i++)
          {
          stopcond+=(zetasold[i]-zetas0[i])*(zetasold[i]-zetas0[i]);
          }
       stopcond=sqrt(stopcond/numfiles);

       std::cout << stopcond << "  (speedup=" << speedup <<")" << std::endl;

       if(stopcond<minvaluespeedup || fabs(stopcond-stopcondold)<minvaluespeedup) 
         {
         REAL aux=speedup/2;
         speedup=(long int) aux;

         if(speedup<1) 
           {
           speedup=1;
           }
         }

       if(speedup>1 && stopcond<=minvalue)
         {
         stopcond=1.0;
         }      
  
       if(count%4==0)
         {
         for(i=0; i<numfiles; i++)
            {
            temp1=(zetas0[i]-zetasold[i]);
            zetasold[i]+=2.0*temp1;
            }
         count=0;
         std::cout << "overrelaxed" << std::endl;
         }
       else
         {        
         for(i=0; i<numfiles; i++)
            {
            zetasold[i]=zetas0[i];
            }
         }
       count++;
       }

  if(fileexists==1)
    {
    outfile.open("zetas.dat", std::ios::out);
    outfile.precision(18);    // <--- print precision
    for(i=0; i<numfiles; i++)
       {
       outfile << zetasold[i] << " ";
       }
    outfile.close();
    }

  delete [] zetasold;
  }



// compute the partition functions of the bootstrapped data
void Data::computezetas(void) 
  {
  const REAL minvalue=1.0e-6;
  const REAL minvaluespeedup=1.0e-4;

  REAL stopcond, stopcondold, energy, temp1, temp2;
  int i, k, j, l, first;
  REAL *zetasold;
  int speedup;
  long int count=1;

  
  std::cout << "\nSelfconsistent computation started\n";
  std::cout << "remainder (stopvalue="<<minvalue<<"):" << std::endl;

 
  zetasold = new REAL [numfiles];
  for(i=0; i<numfiles; i++)
     {
     zetasold[i]=zetas0[i];
     }

  speedup=1;
  stopcond=1.0;
  while(stopcond>minvalue)
       {
       for(k=0; k<numfiles; k++) // estimate zetas[k]
          {
          first=0;

             for(j=0; j<campione_loc; j+=speedup)    //  this means for all the points 
                {                                 

                if(j+speedup<campione_loc)
                  {     
                  energy=0.0;
                  for(l=0; l<speedup; l++)
                     {
                     energy+=7.0*volume*6.0*(1.0-dati[j+l][0]);
                     }
                  energy/=speedup;
                  }
                else
                  {
                  energy=7.0*volume*6.0*(1.0-dati[j][0]);
                  }

                temp1=log(sample_loc[0]/speedup)-zetasold[0]+(betas[k]-betas[0])*energy;
                for(l=1; l<numfiles; l++)
                   {
                   temp2=log(sample_loc[l]/speedup)-zetasold[l]+(betas[k]-betas[l])*energy;
                   temp1=logsum(temp1, temp2);
                   }

                if(first==0)
                  {
                  zetas[k]=-temp1;
                  } 
                else
                  {
                  zetas[k]=logsum(zetas[k], -temp1);
                  }
                
                first=1;
                }
          }
      
       stopcondold=stopcond;
       stopcond=0.0;
       for(i=0; i<numfiles; i++)
          {
          stopcond+=(zetasold[i]-zetas[i])*(zetasold[i]-zetas[i]);
          }
       stopcond=sqrt(stopcond/numfiles);

       std::cout << stopcond << "  (speedup=" << speedup <<")" << std::endl;
 
       if(stopcond<minvaluespeedup || fabs(stopcond-stopcondold)<minvaluespeedup) 
         {
         REAL aux=speedup/2;
         speedup=(long int) aux;
         if(speedup<1) 
           {
           speedup=1;
           }
         }

       if(speedup>1 && stopcond<=minvalue)
         {
         stopcond=1.0;
         }      

       if(count%4==0)
         {
         for(i=0; i<numfiles; i++)
            {
            temp1=(zetas[i]-zetasold[i]);
            zetasold[i]+=2.0*temp1;
            }
         count=0;
         std::cout << "overrelaxed" << std::endl;
         }
       else
         {        
         for(i=0; i<numfiles; i++)
            {
            zetasold[i]=zetas[i];
            }
         }
       count++;

       }

  delete [] zetasold;
  }



void Data::computeobs(REAL beta)
  {
  int first;
  long int f, j, k;
  REAL energy, Z, temp1, temp2;
  REAL *h;
  const int N=2;

  h=new REAL [N*numobs];

  Z=0.0; // just to avoid warning during compilation!

  // sum up all data
  first=0;
  for(f=0; f<campione_loc; f++)              //
     {                                       //  for all data

     energy=7.0*volume*6.0*(1.0-dati[f][0]);
     temp1=log(sample_loc[0])-zetas[0]+(beta-betas[0])*energy;
     for(j=1; j<numfiles; j++)
        {
        temp2=log(sample_loc[j])-zetas[j]+(beta-betas[j])*energy;
        temp1=logsum(temp1, temp2);
        }     

     if(first==0)
       {
       Z=-temp1;
       for(k=0; k<numobs; k++)
          {
          h[N*k]   =     log(dati[f][k])-temp1; // Q
          h[N*k+1] = 2.0*log(dati[f][k])-temp1; // Q^2
          }
       }
     else
       {
       Z=logsum(Z, -temp1);
       for(k=0; k<numobs; k++)
          {
          h[N*k]  =logsum(h[N*k],       log(dati[f][k])-temp1);
          h[N*k+1]=logsum(h[N*k+1], 2.0*log(dati[f][k])-temp1);
          }
       } 
     first=1;

     }

  // normalize helper
  for(j=0; j<N*numobs; j++)
     {
     temp1=h[j]-Z;
     h[j]=exp(temp1);
     }

  //====================
  // compute observables
  //====================
  // first observable: mean plaquette
  medie[0] = h[N*0]-energy_off; // mean 
  medie[1] = volume*(h[N*0+1]-h[N*0]*h[N*0]); // susceptivity
  // second observable: polyakov
  medie[2] = h[N*1]; // mean 
  medie[3] = volume*(h[N*1+1]-h[N*1]*h[N*1]); // susceptibility

  delete [] h;
  }



REAL Data::betapc(REAL beta_low, REAL beta_high, REAL precision)
  {
  int i;
  REAL beta[4];
  REAL suscpoly[4];

  beta[0]=beta_low;
  beta[3]=beta_high;

  if(beta[3]<beta[0])
    {
    return -1.0;
    }

  std::cout.precision(13);    // <--- print precision

  std::cout << "started computation for beta="<<beta[3]<<"...";
  // compute observables
  computeobs(beta[3]);
  suscpoly[3]=medie[3]; // medie[3] is the polyakov loop susceptibility
  std::cout << " terminated (value=" << medie[3] << ")" << std::endl;


  std::cout << "started computation for beta="<<beta[0]<<"...";
  // compute observables
  computeobs(beta[0]);
  suscpoly[0]=medie[3]; // medie[3] is the polyakov loop susceptibility
  std::cout << " terminated (value=" << medie[3] << ")" << std::endl;

  while(beta[3]-beta[0]> precision)
       {
       beta[1]=beta[0]+    (beta[3]-beta[0])/3.0;
       beta[2]=beta[0]+2.0*(beta[3]-beta[0])/3.0;

       for(i=1; i<3; i++)
          {
          std::cout << "started computation for beta="<<beta[i]<<"...";
          // compute observables
          computeobs(beta[i]);
          suscpoly[i]=medie[3]; // medie[3] is the polyakov loop susceptibility
          std::cout << " terminated (value=" << medie[3] << " " << "delta="<< beta[3]-beta[0] << ")" << std::endl;
          }

       if(suscpoly[0]<suscpoly[1] && suscpoly[2]<suscpoly[1]) // maximum near beta[1]
         {
         beta[3]=beta[2];
         suscpoly[3]=suscpoly[2];
         }
       else 
         {
          // maximum near beta[2]
         beta[0]=beta[1];
         suscpoly[0]=suscpoly[1];
         }
       }
  std::cout << "Max. position: beta=" << (beta[3]+beta[0])/2.0 << " +- " << (beta[3]-beta[0])/2.0 << std::endl;

  return (beta[3]+beta[0])/2.0;
  }



/* ----- auxilliary functions ------- */

// compute log(e^a+e^b)
// =a+log(1+e^(b-a)) if a>b 
// =b+log(1+e^(a-b)) if b>a
REAL logsum(REAL a, REAL b)
     {
     REAL diff, ris;
     diff=a-b;
     if(diff >= 0) 
       {
       ris=log1p(exp(-diff));
       ris+=a;
       }
     else
       {
       ris=log1p(exp(diff));
       ris+=b;
       }
     return ris;
     }


// compute log(e^a-e^b)
// =a+log(1-e^(b-a))
REAL logdiff(REAL a, REAL b)
     {
     REAL diff, ris;

     diff=a-b;
     ris=log1p(-exp(-diff));
     ris+=a;
     return ris;
     }


/* ------ MAIN ------ */


int main(int argc, char **argv)
   {
   const int boot_iter=10;
   char *in_file;
   int iter;
   REAL bootsample[boot_iter];
   REAL ris, err;
   std::ofstream file1; 
 

   if(argc != 2)
     {
     std::cout << "Usage: "<< argv[0] << " input_file" <<"\n\n";\
     std::cout << "INPUT FILE TEMPLATE\n";
     std::cout << "outfile  ris.dat\n";
     std::cout << "betamin  0.670\n";
     std::cout << "betamax  0.705\n";
     std::cout << "betastep 0.001\n"; 
     std::cout << "volume   42875\n";
     std::cout << "numfiles 5  #number of data files\n";
     std::cout << "./DATA/dati_35_0.670.dat    0.670    128676    32  #filename  beta  sample  block\n";
     std::cout << "./DATA/dati_35_0.675.dat    0.675    127475    36\n";
     std::cout << "./DATA/dati_35_0.680.dat    0.680    109122    42\n";
     std::cout << "./DATA/dati_35_0.685.dat    0.685    105498    42\n";
     std::cout << "./DATA/dati_35_0.690.dat    0.690    76747     48\n\n";
     std::cout << "IMPORTANT: the first column of the data files has to be the energy density!\n\n";
     std::cout << "The results are <O> and volume*(<O^2>-<O>^2)\n\n";


     return 0;
     }
   else
     {
     in_file=argv[1];
     }


   // initialize parameters
   if(readinput(in_file)!=0)
     {
     return 1;
     }

   // initialize random number generator
   initrand(0);

   // this block is important: 
   // the Data object must be destroied while 
   //  the global variable are still defined !!!
   {
   Data x;

   // initialize data
   x.initfromfile(dati_in);

   // compute the partition functions of the original data
   x.computezetas0();

   for(iter=0; iter<boot_iter; iter++)
      {
      std::cout << "\nSTARTED BOOTSTRAP NUMBER " << iter <<std::endl; 

      // generate bootstrapped data
      x.bootstrapdata();

      // compute the partition functions of the bootstrapped data
      x.computezetas();
     
      bootsample[iter]=x.betapc(betamin, betamax, 1.0e-6);

      file1.open (dati_out.c_str(), std::ios::app);
      file1.precision(13);    // <--- print precision
      file1 << iter << "  " << bootsample[iter] << "\n";
      file1.close();
      }

   ris=0.0;
   for(iter=0; iter<boot_iter; iter++)
      {
      ris+=bootsample[iter];
      }
   ris/=boot_iter;

   err=0.0;
   for(iter=0; iter<boot_iter; iter++)
      {
      err+=(ris-bootsample[iter])*(ris-bootsample[iter]);
      }
   err/=boot_iter;
   err=sqrt(err);

   /*
   std::cout << "FINAL RESULT: beta_pc = " << ris << " +/- " << err << std::endl;

   file1.open (dati_out.c_str(), std::ios::app);
   file1.precision(13);    // <--- print precision
   file1 << "## FINAL:  " << ris << "  " << err;
   file1.close();
   */

   }

   // delete the global variables
   delete [] dati_in;
   delete [] block;
   delete [] sample;
   delete [] numblocks;
   delete [] betas;

   return 0;
   }
