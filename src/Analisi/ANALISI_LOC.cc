#include<cmath>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include<new>
#include<string>

std::string dati_in;
std::string dati_out;
long int block;
long int sample;
int numobs;
long int numblocks;


double myround(double r)
    {
    if(r>0)
      {
      return floor(r + 0.5);
      }
    else
      {
      return ceil(r - 0.5);
      }
    }


int readinput(char *in_file)
    {
    std::ifstream input;
    std::string str, str2;
    long int temp_li;
    const int char_to_ignore=500;

    input.open(in_file, std::ios::in);
    if(input.is_open())
      {
      while(input.good())
           {
           input >> str;

           if(str=="datafile")
             { 
             input >> str2;
             dati_in=str2;
             }
           else if(str=="outfile")
                  {
                  input >> str2;
                  dati_out=str2;
                  }
           else if(str=="block")
                  {
                  input >> temp_li;
                  block=temp_li;
                  }
           else if(str=="sample") 
                  {
                  input >> temp_li;
                  sample=temp_li;
                  }
           else if(str=="numobs") 
                  {
                  input >> temp_li;
                  numobs=temp_li;
                  }
           else
               {
               std::cerr << "Error: unrecognized option \"" << str << "\" in the file \"" << in_file << "\"\n";
               return 1;     
               }

           input.ignore(char_to_ignore, '\n');
           }

      input.close();

      // number of blocks
      numblocks=sample/block;

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
  double** dati;
  double*  medie;
  double*  errori;
  double** jacksample;
  int d_latox, d_latoy, d_latoz, d_latot;
  double d_beta;
public:
  Data();
  ~Data();
  void initfromfile(std::string nome_file);
  void jacksamplegen(void);
  void computeobs(void);
  void printobs(std::string nome_file);
}; 


Data::Data(void)
  {
  long int i;
  
  dati = new double * [sample];
  for(i=0; i<sample; i++)
     {
     dati[i] = new double [numobs];
     }

  medie  = new double [3*numobs];
  errori = new double [3*numobs];

  jacksample = new double * [numblocks];
  for(i=0; i<numblocks; i++)
     {
     jacksample[i]=new double [3*numobs];
     }
  }


Data::~Data(void)
  {
  long int i;
  
  for(i=0; i<sample; i++)
     {
     delete [] dati [i];
     }
  delete [] dati;

  delete [] medie;
  delete [] errori;

  for(i=0; i<numblocks; i++)
     {
     delete [] jacksample[i];
     }
  delete [] jacksample;
  }


void Data::initfromfile(std::string nome_file)
  {
  std::ifstream filein;
  int lx, ly, lz, lt;
  double b, temp1, temp2;
  long int i;

  filein.open(nome_file.c_str(), std::ios::in);
  if(filein.is_open())
    {
    filein >> lx;
    filein >> ly;
    filein >> lz;
    filein >> lt;
    filein >> b;

    d_latox=lx;
    d_latoy=ly;
    d_latoz=lz;
    d_latot=lt;
    d_beta=b;

    for(i=0; i<sample; i++)
       {
       filein >> temp1;
       filein >> temp2;
       dati[i][0]=(temp1+temp2)*0.5; // average plaquette

       filein >> temp1;
       filein >> temp2;
       dati[i][1]=sqrt(temp1*temp1+temp2*temp2); // polyakov loop

       filein.ignore(1000, '\n');
       }
    filein.close();
    }
  else
    {
    std::cerr << "Error in opening the file "<<nome_file<<" \n";
    exit(1);
    }
  }


// generate jacknife samples
void Data::jacksamplegen(void)
  {
  long int i, b;
  int j;
  double* h; // for helper
  const double vol=(double) ( d_latox*d_latoy*d_latoz*d_latot );

  h = new double[3*numobs];

  // initialize to zero medie
  for(j=0; j<3*numobs; j++)
     {
     medie[j]=0.0;
     }
  // sum up all data
  for(i=0; i<numblocks*block; i++)
     {
     for(j=0; j<numobs; j++)
        {
        medie[3*j]  +=dati[i][j]; // Q   
        medie[3*j+1]+=dati[i][j]*dati[i][j]; // Q^2
        medie[3*j+2]+=dati[i][j]*dati[i][j]*dati[i][j]*dati[i][j]; // Q^4
        }
     }
  // generate jacknife samples
  for(b=0; b<numblocks; b++)
     {
     // initialize helper
     for(j=0; j<3*numobs; j++)
        {
        h[j]=medie[j];
        }
     // subtract the b-th block
     for(i=b*block; i<(b+1)*block; i++)
        {
        for(j=0; j<numobs; j++)
           {
           h[3*j]  -=dati[i][j]; // Q   
           h[3*j+1]-=dati[i][j]*dati[i][j]; // Q^2
           h[3*j+2]-=dati[i][j]*dati[i][j]*dati[i][j]*dati[i][j]; // Q^4
           } 
        }
     // normalize helper
     for(j=0; j<3*numobs; j++)
        {
        h[j]/=(double) ((numblocks-1)*block);
        }
     // compute observables
     for(j=0; j<numobs; j++)
        {
        jacksample[b][3*j] = h[3*j]; // media
        jacksample[b][3*j+1] = (h[3*j+1]-h[3*j]*h[3*j])*vol; // vol(<Q^2>-<Q>^2)
        jacksample[b][3*j+2] = 1.0-h[3*j+2]/h[3*j+1]/h[3*j+1]/3.0; // 1-<Q^4>/<Q^2>^2/3
        }
     }
  delete [] h;
  }


void Data::computeobs(void)
  {
  long int i;
  int j;
  double temp;

  for(j=0; j<3*numobs; j++)
     {
     temp=0.0;
     for(i=0; i<numblocks; i++)
        {
        temp+=jacksample[i][j];
        }
     temp/=(double)numblocks;
     medie[j]=temp;
     }

  for(j=0; j<3*numobs; j++)
     {
     temp=0.0;
     for(i=0; i<numblocks; i++)
        {
        temp+=(jacksample[i][j]-medie[j])*(jacksample[i][j]-medie[j]);
        }
     temp/=(double)numblocks;
     temp*=(double)(numblocks-1);
     errori[j]=sqrt(temp);
     }
  }


void Data::printobs(std::string nome_file)
  {
  std::ofstream file1; 
  int i;

  file1.open (nome_file.c_str(), std::ios::app);
  file1 << d_beta;

  file1.precision(13);    // <--- print precision

  for(i=0; i<numobs; i++)
     {
     file1 << " " << medie[3*i]   << " " << errori[3*i];
     file1 << " " << medie[3*i+1] << " " << errori[3*i+1];
     file1 << " " << medie[3*i+2] << " " << errori[3*i+2];
     }
  file1 << "\n";
  file1.close();
  }



int main(int argc, char **argv)
   {
   char *in_file;


  if(argc != 2)
    {
    std::cout << "Usage: "<< argv[0] << " input_file" <<"\n";
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

   Data x;

   // initialize data
   x.initfromfile(dati_in);
  
   // generate jacknife samples
   x.jacksamplegen();

   // compute observables
   x.computeobs();

   // print results
   x.printobs(dati_out);

   return 0;
   }
