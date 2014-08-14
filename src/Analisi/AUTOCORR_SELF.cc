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
#include<fstream>
#include<iostream>
#include<new>
#include<string>

// For details see e.g app. C of Madras, Sokal ``The pivot algorithm: A highly efficient Monte Carlo method for the self-avoiding walk''
//                            Journal of statistical physics, vol. 50, pag. 109 (1988)
//                 or U. Wolff ``Monte Carlo errors with less errors'' arXiv:hep-lat/0306017
//	                       Comput.Phys.Commun.156:143-153,2004; Erratum-ibid.176:383,2007

std::string dati_in;
long int sample;
int columns;


class Data {
private:
  double** dati;
  double* medie;
  double** C;
  double* tau;
  long int* block;
public:
  Data();
  ~Data();
  void initfromfile(std::string nome_file);
  void computemean();
  void computeC(int col);
  void computetau(int col);
  void computetau(void);
  void printdataonscreen();
  void printmeanonscreen();
  void printtauonscreen();

}; 


Data::Data()
  {
  long int i;
  dati = new double * [sample];
  for(i=0; i<sample; i++)
     {
     dati[i] = new double [columns];
     }

  medie = new double [columns];
  tau = new double [columns];
  block=new long int [columns];
 
  C = new double * [sample];
  for(i=0; i<sample; i++)
     {
     C[i] = new double [columns];
     }
  }


Data::~Data(void)
  {
  long int i;
  
  for(i=0; i<sample; i++)
     {
     delete [] dati[i];
     }
  delete [] dati;

  delete [] medie;
  delete [] tau;
  delete [] block;

  for(i=0; i<sample; i++)
     {
     delete [] C[i];
     }
  delete [] C;
  }


void Data::initfromfile(std::string nome_file)
  {
  int j;
  std::ifstream filein;
  double temp;
  long int i;

  filein.open(nome_file.c_str(), std::ios::in);
  if(filein.is_open())
    {
    filein >> j;
    filein >> j;
    filein >> j;
    filein >> j;
    filein >> temp;

    for(i=0; i<sample; i++)
       {
       for(j=0; j<columns; j++)
          {
          filein>> temp;
          dati[i][j]=temp;
          }       
       }

    filein.close();
    }
  else
    {
    std::cerr << "Error in opening the file "<<nome_file<<" \n";
    exit(1);
    }
  }


void Data::computemean(void)
  {
  long int i;
  int col;

  for(col=0; col<columns; col++)
     {
     medie[col]=0.0;
     } 


  for(i=0; i<sample; i++)
     {
     for(col=0; col<columns; col++)
        {
        medie[col]+=dati[i][col];
        } 
     }

  for(col=0; col<columns; col++)
     {
     medie[col]/=(double) sample;
     } 
  }


void Data::computeC(int col)
  {
  long int t, i;
  double aux;

  for(t=0; t<block[col]; t++)
     {
     aux=0.0;
     for(i=0; i<sample-block[col]; i++)
        {
        aux+=(dati[i][col]-medie[col])*(dati[i+t][col]-medie[col]);
        }
     aux/=(double)(sample-block[col]);
     C[t][col]=aux;
     }   
  }


void Data::computetau(int col)
  {
  long int t;
  double aux;

  aux=0.0;
  for(t=0; t<block[col]; t++)
     {
     if(t==0)
       {
       aux+=1.0;
       }
     else
       {
       aux+=2.0*C[t][col]/C[0][col];
       }
     } 
  tau[col]=aux;
  }


void Data::computetau(void)
  {
  int col;
  const double S=1.5; // see Wolff paper
  double tau_h;
  double deriv;

  for(col=0; col<columns; col++)
     {
     block[col]=1;
     deriv=1.0;

     while(block[col]<sample/2 && deriv>0.0)
          { 
          block[col]+=1;
          computeC(col);
          computetau(col);

          if(tau[col]<0.5)
            {
            tau_h=0.6;
            }
          else
            {
            tau_h=S/log((2.0*tau[col]+1.0)/(2.0*tau[col]-1.0));
            }
   
          deriv=exp(-((double) block[col])/tau_h)-tau_h/sqrt((double) (block[col]*sample));
          }
     }
  }


void Data::printdataonscreen(void)
  {
  long int i;
  int j;

  for(i=0; i<sample; i++)
     {
     for(j=0; j<columns; j++)
        {
        std::cout << dati[i][j]<< " ";
        }
     std::cout << "\n";
     }
  }


void Data::printmeanonscreen(void)
  {
  int j;

  for(j=0; j<columns; j++)
     {
     std::cout << medie[j]<< " ";
     }
  std::cout << "\n";
  }


void Data::printtauonscreen(void)
  {
  int col;

  for(col=0; col<columns; col++)
     {
     std::cout << "=================\n"; 
     std::cout << "COLUMN "<< col<<"\n";
     if(block[col]<sample/2)
       {
       std::cout << "estimated optimal block   : " << block[col] <<"\n";
       }
     else
       {
       std::cout << "estimated optimal block   : " << block[col] <<" (WARNING: sample may be too small!)\n";
       }
     std::cout << "estimated autocorrelation : " << tau[col] << "    ";
     std::cout << "(" << tau[col]*sqrt(2.0/((double)sample)*(2.0*block[col]+1.0))<<")\n";
     }
  std::cout << "\n";
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
           else if(str=="sample") 
                  {
                  input >> temp_li;
                  sample=temp_li;
                  }
           else if(str=="columns") 
                  {
                  input >> temp_li;
                  columns=temp_li;
                  }
           else
               {
               std::cerr << "Error: unrecognized option \"" << str << "\" in the file \"" << in_file << "\"\n";
               return 1;     
               }

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

   // initializa data
   x.initfromfile(dati_in);

   // compute the averages
   x.computemean();  

   // compute autocorrelation times
   x.computetau();

   // print autocorrelation times on screen
   x.printtauonscreen();

   return 0;
   }
