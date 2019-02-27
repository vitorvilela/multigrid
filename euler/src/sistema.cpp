//SISTEMA

#include <cstdlib>
#include <iostream>

using std::cout;

#include <stdio.h>
#include <math.h>

#include "sistema.h"
#include "solvers.h"
#include "utilitarios.h"
#include "exatas.h"


void ME1(double ****Ro, int *nx, int *ny, int *nz, double *h, int K, char imp, int niveis)
{
     int t, i, j, k;
     double x, y, z;
     
     int f = 1;
     double e = 2;
     
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
 	      sprintf(titulo_ro, "Ro1_%03i.dat", K);
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



void ME2(double ****Ro, int *nx, int *ny, int *nz, double *h, int K, char imp, int niveis)
{
     int t, i, j, k;
     double a, b, c;
     const double PI = atan(1)*4;
	
     for(t=0; t<niveis; t++)
     {
          for(i=0; i<nx[t]; i++)
          {
               for(j=0; j<ny[t]; j++)
               {
	                for(k=0; k<nz[t]; k++)
				    {
                         if(k < nz[t]/2)
                              Ro[t][i][j][k] = 1 + K;                              
                         else
                              Ro[t][i][j][k] = 1;
			        }
               }
          }
     }
	
	 if(imp == 's')
	 { 
          FILE *sai_ro;
          char titulo_ro[10];
 	      sprintf(titulo_ro, "Ro2_%03i.dat", K);
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


void ME3(double ****Ro, int *nx, int *ny, int *nz, double *h, int K, char imp, int niveis)
{
     int t, i, j, k;
     double a, b, c;
     const double PI = atan(1)*4;
     
     int nxi, nxf, nyi, nyf;
     double xi, xf, yi, yf;
     xi = yi = 0.4;
     xf = yf = 0.6;
     	
     for(t=0; t<niveis; t++)
     {          
          nxi = nyi = ceil(xi/h[t]);
          nxf = nyf = ceil(xf/h[t]);
              
          for(i=0; i<nx[t]; i++)          
               for(j=0; j<ny[t]; j++)
                    for(k=0; k<nz[t]; k++)	                             
                         Ro[t][i][j][k] = 1;
                         
          for(i=nxi; i<=nxf; i++)
               for(j=nyi; j<=nyf; j++)
                    for(k=0; k<nz[t]; k++)
                         Ro[t][i][j][k] = 1 + K;
     }
	
	 if(imp == 's')
	 { 
          FILE *sai_ro;
          char titulo_ro[10];
 	      sprintf(titulo_ro, "Ro3_%03i.dat", K);
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


void ME4(double ****Ro, int *nx, int *ny, int *nz, double *h, int K, char imp, int niveis)
{
     int t, i, j, k;
     double a, b, c;
     const double PI = atan(1)*4;
     
     int nxi, nxf, nyi, nyf, nzi, nzf;
     double xi, xf, yi, yf, zi, zf;
     xi = yi = zi = 0.4;
     xf = yf = zf = 0.6;
     	
     for(t=0; t<niveis; t++)
     {          
          nxi = nyi = nzi = ceil(xi/h[t]);
          nxf = nyf = nzf = ceil(xf/h[t]);
              
          for(i=0; i<nx[t]; i++)          
               for(j=0; j<ny[t]; j++)
                    for(k=0; k<nz[t]; k++)	                             
                         Ro[t][i][j][k] = 1;
                         
          for(i=nxi; i<=nxf; i++)
               for(j=nyi; j<=nyf; j++)
                    for(k=nzi; k<nzf; k++)
                         Ro[t][i][j][k] = 1 + K;
     }
	
	 if(imp == 's')
	 { 
          FILE *sai_ro;
          char titulo_ro[10];
 	      sprintf(titulo_ro, "Ro4_%03i.dat", K);
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


void precondicionador(double ****Mp, double ****Mw, double ****Me, double ****Mn, double ****Ms, double ****Mu, double ****Md, double ****Ro, double *h, int *nx, int *ny, int *nz, int niveis, char imp)
{
     int t, i, j, k;
     
     double A;
     
     for(t=0; t<niveis; t++)
          for(i=1; i<(nx[t]-1); i++)
               for(j=1; j<(ny[t]-1); j++)
			        for(k=1; k<(nz[t]-1); k++)
                    {       
                         A = -(1/Ro[t][i][j][k])*1/(h[t]*h[t]);
                         
                         Mn[t][i][j][k] = -(0.5*(1/(Ro[t][i][j+1][k])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         if(Mn[t][i][j][k] != A) 
                              Mn[t][i][j][k] = 0;                         
                                                 
                         Ms[t][i][j][k] = -(0.5*(1/(Ro[t][i][j-1][k])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         if(Ms[t][i][j][k] != A) 
                              Ms[t][i][j][k] = 0;                         
                         
                         Mw[t][i][j][k] = -(0.5*(1/(Ro[t][i-1][j][k])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         if(Mw[t][i][j][k] != A) 
                              Mw[t][i][j][k] = 0;                         
                                                  
                         Me[t][i][j][k] = -(0.5*(1/(Ro[t][i+1][j][k])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         if(Me[t][i][j][k] != A) 
                              Me[t][i][j][k] = 0;
                         
                         Mu[t][i][j][k] = -(0.5*(1/(Ro[t][i][j][k+1])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         if(Mu[t][i][j][k] != A) 
                              Mu[t][i][j][k] = 0;
                         
                         Md[t][i][j][k] = -(0.5*(1/(Ro[t][i][j][k-1])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         if(Md[t][i][j][k] != A) 
                              Md[t][i][j][k] = 0;
                         
                         Mp[t][i][j][k] = -(Me[t][i][j][k]+Mw[t][i][j][k]+Mn[t][i][j][k]+Ms[t][i][j][k]+Mu[t][i][j][k]+Md[t][i][j][k]);
                    } 
                    
     if(imp == 's')
	 { 
          FILE *sai_c;
          char titulo_c[10];
 	      sprintf(titulo_c, "Coef.dat");
       	  sai_c = fopen(titulo_c, "w");
          fprintf(sai_c,"variables = \"x\" , \"y\" , \"z\" , \"C\" \n");
          fprintf(sai_c,"zone i=%04d j=%04d k=%04d f=point \n", nx[niveis-1]-2, ny[niveis-1]-2, nz[niveis-1]-2);
        	
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
                         fprintf(sai_c, "%e   %e   %e   %e\n", vadx, vady, vadz, Mw[niveis-1][i][j][k]);
                         vadz += h[niveis-1];
				    }
				    vadx += h[niveis-1];
                }
                vady += h[niveis-1];
          }
          fclose(sai_c);
     }   
}



     
