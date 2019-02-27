//funcoes.cpp

#include <cstdlib>
#include <iostream>

using std::cout;

#include <stdio.h>
#include <math.h>

#include "sistema.h"
#include "solvers.h"
#include "utilitarios.h"
#include "funcoes.h"



void fronteira(double ***P, double h, int nx, int ny, int nz, double D, double L, double S, int fun, char fron, int sit)
{
     if(fun == 1)
     {
          if(fron == 'd')
               D1(P, h, nx, ny, nz, D, L, S, sit);
          else if(fron == 'n')          
               N1(P, h, nx, ny, nz, D, L, S, sit);                 
     }
}

void E1(double ***Pex, double h, int nx, int ny, int nz, char imp)
{         
     int nro = 1;
     int i, j, k;
     double x, y, z;
     
     double Pi = 4*atan(1);
     
     for(i=0; i<nx; i++)
     {
          x = 0.5*h + (i-1)*h;
          for(j=0; j<ny; j++)
          {
               y = 0.5*h + (j-1)*h;			
		       for(k=0; k<nz; k++)
		       {
                    z = 0.5*h + (k-1)*h;
				    Pex[i][j][k] = 2*x*x - 5*y*y + 4*z*z -1/3;
		       }
          }
     }
     
     if(imp == 's')
     {
          FILE *sai_e;
          char titulo_e[10];
          sprintf(titulo_e, "Pex_%03i.dat", nro);
       	  sai_e = fopen(titulo_e, "w");
          fprintf(sai_e,"variables = \"x\" , \"y\" , \"z\" , \"Pex\" \n");
          fprintf(sai_e,"zone i=%04d j=%04d k=%04d f=point \n", nx-2, ny-2, nz-2);
        	
		  double vadx, vady, vadz;
          
          vady = 0.0;
          for(j=1; j<(ny-1); j++)
          {
               vadx = 0.0;
               for(i=1; i<(nx-1); i++)
               {
                    vadz = 0.0;
				    for(k=1; k<(nz-1); k++)
				    {
                         fprintf(sai_e, "%e   %e   %e   %e\n", vadx, vady, vadz, Pex[i][j][k]);
                         vadz += h;
				    }
			        vadx += h;
               }
               vady += h;
          }
          fclose(sai_e);
     }     
}

void D1(double ***P, double h, int nx, int ny, int nz, double D, double L, double S, int sit)
{
     int i, j, k;
     double x, y, z;
     
     for(j=0; j<ny; j++)
     {
          y = 0.5*h + (j-1)*h;              
	      for(k=0; k<nz; k++)
	      {
	           z = 0.5*h + (k-1)*h;
               if(sit == 1)     
               {               
                    P[0][j][k] = 2*(-5*y*y+4*z*z-1/3) - P[1][j][k];
                    P[nx-1][j][k] = 2*(2*L*L-5*y*y+4*z*z-1/3) - P[nx-2][j][k];
		       }
               else    
        	   {
                    P[0][j][k] = -P[1][j][k];
                    P[nx-1][j][k] = -P[nx-2][j][k];
               }		
          }
     }

     for(i=0; i<nx; i++)
     {
          x = 0.5*h + (i-1)*h;
	      for(k=0; k<nz; k++)
	      {
		       z = 0.5*h + (k-1)*h;              
        	   if(sit == 1)
               {
                    P[i][0][k] = 2*(2*x*x+4*z*z-1/3) - P[i][1][k];
                    P[i][ny-1][k] = 2*(2*x*x+4*z*z-5*D*D-1/3) - P[i][ny-2][k];
               }
               else 
               {
                    P[i][0][k] = -P[i][1][k];
                    P[i][ny-1][k] = -P[i][ny-2][k];
           	   }
	      }
     } 
     
     for(i=0; i<nx; i++)
     {
          x = 0.5*h + (i-1)*h;
	      for(j=0; j<ny; j++)
	      {
               y = 0.5*h + (j-1)*h;
               if(sit == 1) 
               {
                    P[i][j][0] = 2*(2*x*x-5*y*y-1/3) - P[i][j][1];
                    P[i][j][nz-1] = 2*(2*x*x-5*y*y+4*S*S-1/3) - P[i][j][nz-2];
               }
               else 
               {
              		P[i][j][0] = -P[i][j][1];
               		P[i][j][nz-1] = -P[i][j][nz-2];
               }
	      }
     } 
}

void N1(double ***P, double h, int nx, int ny, int nz, double D, double L, double S, int sit)
{
     int i, j, k;
     double x, y, z;
     
     if(sit == 1)
          normaliza(P, nx, ny, nz, h);
     
     for(j=0; j<ny; j++)
     {
          y = 0.5*h + (j-1)*h;              
	      for(k=0; k<nz; k++)
	      {
	           z = 0.5*h + (k-1)*h;
               if(sit == 1)     
               {               
                    P[0][j][k] = P[1][j][k]; 
                    P[nx-1][j][k] = P[nx-2][j][k] + h*(4*L);
		       }
               else    
        	   {
                    P[0][j][k] = P[1][j][k];
                    P[nx-1][j][k] = P[nx-2][j][k];
               }		
          }
     }

     for(i=0; i<nx; i++)
     {
          x = 0.5*h + (i-1)*h;
	      for(k=0; k<nz; k++)
	      {
		       z = 0.5*h + (k-1)*h;              
        	   if(sit == 1)
               {
                    P[i][0][k] = P[i][1][k];
                    P[i][ny-1][k] = P[i][ny-2][k] + h*(-10*D);
               }
               else 
               {
                    P[i][0][k] = P[i][1][k];
                    P[i][ny-1][k] = P[i][ny-2][k];
           	   }
	      }
     } 
     
     for(i=0; i<nx; i++)
     {
          x = 0.5*h + (i-1)*h;
	      for(j=0; j<ny; j++)
	      {
               y = 0.5*h + (j-1)*h;
               if(sit == 1) 
               {
                    P[i][j][0] = P[i][j][1];
                    P[i][j][nz-1] = P[i][j][nz-2] + h*(8*S);
               }
               else 
               {
              		P[i][j][0] = P[i][j][1];
               		P[i][j][nz-1] = P[i][j][nz-2];
               }
	      }
     } 
}



