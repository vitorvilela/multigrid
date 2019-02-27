//sistema.cpp

#include <cstdlib>
#include <iostream>
using std::cout;

#include <stdio.h>
#include <math.h>


void ME(double ****Ro, int *nx, int *ny, int *nz, double *h, int K, char imp, int niveis)
{
     int t, i, j, k;
     double x, y, z;
     
     int f = 1;
     double e = 20;
     
     const double Pi = atan(1)*4;
	
     for(t=0; t<niveis; t++)
     {
          for(i=0; i<nx[t]; i++)
          {
               x = 0.5*h[t] + (i-1)*h[t];
               
               for(j=0; j<ny[t]; j++)
               {
                    y = 0.5*h[t] + (j-1)*h[t];    
                        
	                for(k=0; k<nz[t]; k++)
				    {
                         z = 0.5*h[t] + (k-1)*h[t];    
                            					     
                         Ro[t][i][j][k] = 1 + K*pow(sin(f*Pi*x)*sin(f*Pi*y)*sin(f*Pi*z),e);                         
			        }
               }
          }
     }
	
	 if(imp == 's')
	 { 
          FILE *sai_ro;
          char titulo_ro[10];
 	      sprintf(titulo_ro, "Ro_%05i.dat", K);
       	  sai_ro = fopen(titulo_ro, "w");
          fprintf(sai_ro,"variables = \"x\" , \"y\" , \"z\" , \"Ro\" \n");
          fprintf(sai_ro,"zone i=%04d j=%04d k=%04d f=point \n", nx[niveis-1]-2, ny[niveis-1]-2, nz[niveis-1]-2);
        	
		  double vadx, vady, vadz;
		  
          vady = 0.0;
          for(j=1; j<(ny[niveis-1]-1); j++)
          {
               vadx = 0.0;
               for(i=1; i<(nx[niveis-1]-1); i++)
               {
		            vadz = 0.0;
				    for(k=1; k<(nz[niveis-1]-1); k++)
				    {
                         fprintf(sai_ro, "%e   %e   %e   %e\n", vadx, vady, vadz, Ro[niveis-1][i][j][k]);
                         vadz += h[niveis-1];
				    }
				    vadx += h[niveis-1];
                }
                vady += h[niveis-1];
          }
          fclose(sai_ro);
     }   
}


void coeficientes(double ****Ap, double ****Aw, double ****Ae, double ****An, double ****As, double ****Au, double ****Ad, double ****Ro, double *h, int *nx, int *ny, int *nz, int niveis)
{
     int t, i, j, k;
     
     for(t=0; t<niveis; t++)
          for(i=1; i<(nx[t]-1); i++)
               for(j=1; j<(ny[t]-1); j++)
			        for(k=1; k<(nz[t]-1); k++)
                    {       
                         An[t][i][j][k] = -(0.5*(1/(Ro[t][i][j+1][k])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         As[t][i][j][k] = -(0.5*(1/(Ro[t][i][j-1][k])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         Aw[t][i][j][k] = -(0.5*(1/(Ro[t][i-1][j][k])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         Ae[t][i][j][k] = -(0.5*(1/(Ro[t][i+1][j][k])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         Au[t][i][j][k] = -(0.5*(1/(Ro[t][i][j][k+1])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         Ad[t][i][j][k] = -(0.5*(1/(Ro[t][i][j][k-1])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         Ap[t][i][j][k] = -(Ae[t][i][j][k]+Aw[t][i][j][k]+An[t][i][j][k]+As[t][i][j][k]+Au[t][i][j][k]+Ad[t][i][j][k]);
                    } 
}


void f(double ***B, double ***Pex, double ***Ap, double ***Aw, double ***Ae, double ***An, double ***As, double ***Au, double ***Ad, int nx, int ny, int nz)
{
     int t, i, j, k;
              
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               for(k=1; k<(nz-1); k++)	
                    B[i][j][k] = Ap[i][j][k]*Pex[i][j][k]+Ae[i][j][k]*Pex[i+1][j][k]+
                                 Aw[i][j][k]*Pex[i-1][j][k]+An[i][j][k]*Pex[i][j+1][k]+
                                 As[i][j][k]*Pex[i][j-1][k]+Au[i][j][k]*Pex[i][j][k+1]+
					             Ad[i][j][k]*Pex[i][j][k-1];	
}
     
                 

     
