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

#define REAL long double


////// global variables 
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

//////

const int numobs=2;
const int numobs_tot_jack=4;

const REAL pi=3.141592653589793238462643383279502884197169399375105820974944;

const REAL energy_off=0.0;

//////



REAL logsum(REAL a, REAL b);
REAL logdiff(REAL a, REAL b);





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
  REAL*** dati;
  REAL*  medie;
  REAL*  errori;
  REAL** jacksample;
  REAL* zetas;
public:
  Data();
  ~Data();
  void initfromfile(std::string *nome_file);
  void computezetas(void);
  void jacksamplegen(REAL beta);
  void computeobs(void);
  void printobs(std::string nome_file, REAL beta);
}; 



Data::Data(void)
  {
  long int i, j;
  
  dati = new REAL ** [numfiles];
  for(i=0; i<numfiles; i++)
     {
     dati[i] = new REAL * [sample[i]];
     for(j=0; j<sample[i]; j++)
        {
        dati[i][j]=new REAL [numobs];
        }
     }

  medie  = new REAL [numobs_tot_jack];
  errori = new REAL [numobs_tot_jack];

  jacksample = new REAL * [numblockstot];
  for(i=0; i<numblockstot; i++)
     {
     jacksample[i]=new REAL [numobs_tot_jack];
     }

  zetas = new REAL [numfiles];
  }


Data::~Data(void)
  {
  long int i, j;
  
  for(i=0; i<numfiles; i++)
     {
     for(j=0; j<sample[i]; j++)
        {
        delete [] dati[i][j];
        }
     delete [] dati[i];
     }
  delete [] dati;

  delete [] medie;
  delete [] errori;

  for(i=0; i<numblockstot; i++)
     {
     delete [] jacksample[i];
     }
  delete [] jacksample;

  delete [] zetas;
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
               dati[f][i][j]=temp1 + energy_off; // so that it is >0
               }
             else
               {
               dati[f][i][j]=temp1;
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



// compute the partition functions
void Data::computezetas(void) 
  {
  const REAL minvalue=1.0e-7;
  const REAL minvaluespeedup=1.0e-4;

  REAL stopcond, stopcondold, energy, temp1, temp2;
  int i, k, f, j, l, first;
  REAL *zetasold;
  int speedup;
  int count=1;
  
  int fileexists;
  std::ifstream infile;
  std::ofstream outfile;

  std::cout << "\nSelfconsistent computation started\n";
  std::cout << "remainder (stopvalue="<<minvalue<<"):\n";

  zetasold = new REAL [numfiles];

  srand((unsigned)time(0)); 

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
       zetasold[i]=1.0+rand()/RAND_MAX;
       }
    speedup=500;
    fileexists=1;
    }

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
                     energy+=volume*6.0*(1.0-dati[f][j+l][0]);
                     }
                  energy/=speedup;
                  }
                else
                  {
                  energy=volume*6.0*(1.0-dati[f][j][0]);
                  }
                temp1=log(sample[0]/speedup)-zetasold[0]+(betas[k]-betas[0])*energy;
                for(l=1; l<numfiles; l++)
                   {
                   temp2=log(sample[l]/speedup)-zetasold[l]+(betas[k]-betas[l])*energy;
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




// generate jacknife samples
void Data::jacksamplegen(REAL beta)
  {
  int first;
  long int f, i, j, k, b, btot;
  REAL energy, Z, temp1, temp2;
  REAL *h, *h0, Zh;
  const int N=2;

  h0=new REAL [N*numobs];
  h=new REAL [N*numobs];

  Z=0.0; // just to avoid warning during compilation!

  // sum up all data
  first=0;
  for(f=0; f<numfiles; f++)                  //
     {                                       //  for all data
     for(i=0; i<numblocks[f]*block[f]; i++)  //  
        {                                    //

        energy=volume*6.0*(1.0-dati[f][i][0]);
        temp1=log(sample[0])-zetas[0]+(beta-betas[0])*energy;
        for(j=1; j<numfiles; j++)
           {
           temp2=log(sample[j])-zetas[j]+(beta-betas[j])*energy;
           temp1=logsum(temp1, temp2);
           }     

        if(first==0)
          {
          Z=-temp1;
          for(k=0; k<numobs; k++)
             {
             h0[N*k]   =     log(dati[f][i][k])-temp1; // Q
             h0[N*k+1] = 2.0*log(dati[f][i][k])-temp1; // Q^2
             }
          }
        else
          {
          Z=logsum(Z, -temp1);
          for(k=0; k<numobs; k++)
             {
             h0[N*k]  =logsum(h0[N*k],       log(dati[f][i][k])-temp1);
             h0[N*k+1]=logsum(h0[N*k+1], 2.0*log(dati[f][i][k])-temp1);
             }
          } 
        first=1;

        }
     }

 
  // generate jacknife samples         
  btot=0;                              //
  for(f=0; f<numfiles; f++)            // for all blocks
     {                                 //
     for(b=0; b<numblocks[f]; b++)     //
        {

        // initialize helpers
        Zh=Z;
        for(j=0; j<N*numobs; j++)
           {
           h[j]=h0[j];
           }

        // subtract the b-th block
        for(i=b*block[f]; i<(b+1)*block[f]; i++)
           {
           energy=volume*6.0*(1.0-dati[f][i][0]);
           temp1=log(sample[0])-zetas[0]+(beta-betas[0])*energy;
           for(j=1; j<numfiles; j++)
              {
              temp2=log(sample[j])-zetas[j]+(beta-betas[j])*energy;
              temp1=logsum(temp1, temp2);
              }     
 
           Zh=logdiff(Zh, -temp1);
           for(k=0; k<numobs; k++)
             {
             h[N*k]  =logdiff(h[N*k],       log(dati[f][i][k])-temp1);
             h[N*k+1]=logdiff(h[N*k+1], 2.0*log(dati[f][i][k])-temp1);
             }
           }

        // normalize helper
        for(j=0; j<N*numobs; j++)
           {
           temp1=h[j]-Zh;
           h[j]=exp(temp1);
           }

        //====================
        // compute observables
        //====================
        // first observable: mean plaquette
        jacksample[btot][0] = h[N*0]-energy_off; // mean 
        jacksample[btot][1] = volume*(h[N*0+1]-h[N*0]*h[N*0]); // susceptivity
        // second observable: polyakov
        jacksample[btot][2] = h[N*1]; // mean 
        jacksample[btot][3] = volume*(h[N*1+1]-h[N*1]*h[N*1]); // susceptibility

        btot++;       

        }
     }

  delete [] h;
  delete [] h0;
  }


void Data::computeobs(void)
  {
  long int i;
  int j;
  double temp;

  for(j=0; j<numobs_tot_jack; j++)
     {
     temp=0.0;
     for(i=0; i<numblockstot; i++)
        {
        temp+=jacksample[i][j];
        }
     temp/=(double)numblockstot;
     medie[j]=temp;
     }

  for(j=0; j<numobs_tot_jack; j++)
     {
     temp=0.0;
     for(i=0; i<numblockstot; i++)
        {
        temp+=(jacksample[i][j]-medie[j])*(jacksample[i][j]-medie[j]);
        }
     temp/=(double)numblockstot;
     temp*=(double)(numblockstot-1);
     errori[j]=sqrt(temp);
     }
  }


void Data::printobs(std::string nome_file, REAL beta)
  {
  std::ofstream file1; 
  int i;

  file1.open (nome_file.c_str(), std::ios::app);

  file1.precision(13);    // <--- print precision
  file1<< beta;
  for(i=0; i<numobs_tot_jack; i++)
     {
     file1 << " " << medie[i]  << " " << errori[i];
     }
  
  file1<<"\n";

  file1.close();
  }



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



int main(int argc, char **argv)
   {
   char *in_file;
   REAL beta;

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

   // this block is important: 
   // the Data object must be destroied while 
   //  the global variable are still defined !!!
   {
   Data x;

   // initialize data
   x.initfromfile(dati_in);

   // compute the partition functions
   x.computezetas();

   for(beta=betamin; beta<betamax; beta+=betastep)
      {
      std::cout << "started computation for beta="<<beta<<"...";
      // generate jacknife samples
      x.jacksamplegen(beta);

      // compute observables
      x.computeobs();

      // print results
      x.printobs(dati_out, beta);
      std::cout << " terminated\n";
      }

   }

   // delete the global variables
   delete [] dati_in;
   delete [] block;
   delete [] sample;
   delete [] numblocks;
   delete [] betas;

   return 0;
   }
