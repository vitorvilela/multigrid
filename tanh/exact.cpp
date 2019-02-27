//exact.cpp

#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <iostream>
using std::cout;
using std::cin;

#include "system.h"
#include "solvers.h"
#include "tools.h"
#include "exact.h"


void EXACT(double ****Pex, double *h, int *nx, int *ny, int *nz, int levels, char imp)
{         
     int t, i, j, k;
     
     double x, y, z;
     
     double Pi = 4*atan(1);
     
     for(t=0; t<levels; t++)
          for(i=0; i<nx[t]; i++)
          {
               x = 0.5*h[t] + (i-1)*h[t];
               
               for(j=0; j<ny[t]; j++)
               {
                    y = 0.5*h[t] + (j-1)*h[t];			
		            
                    for(k=0; k<nz[t]; k++)
		            {                             
                         z = 0.5*h[t] + (k-1)*h[t];
				         
                         Pex[t][i][j][k] = 2*x*x - 5*y*y + 4*z*z -1/3;
		            }
               }
          }
     
     if(imp == 's')
	 { 
          t = levels-1;  
            
          FILE *out;
          char title[15];
 	      sprintf(title, "Ex %03i.dat", nx[t]-2);
       	  out = fopen(title, "w");
          fprintf(out, "variables = \"x\" , \"y\" , \"z\" , \"Ex\" \n");
          fprintf(out, "zone i=%04d j=%04d k=%04d f=point \n", nx[t]-2, ny[t]-2, nz[t]-2);
        		  
          y = 0.0;
          for(j=1; j<(ny[t]-1); j++)
          {
               x = 0.0;
               for(i=1; i<(nx[t]-1); i++)
               {
		            z = 0.0;
				    for(k=1; k<(nz[t]-1); k++)
				    {
                         fprintf(out, "%e   %e   %e   %e\n", x, y, z, Pex[t][i][j][k]);
                         z += h[t];
				    }
				    x += h[t];
                }
                y += h[t];
          }
          fclose(out);
     }   
}


void SET_BOUNDARY(double ***P, double X, double Y, double Z, double h, int nx, int ny, int nz, char boundaryType, int sit)
{
     if(boundaryType == 'd')
          DIRICHLET(P, X, Y, Z, h, nx, ny, nz, sit);
     else if(boundaryType == 'n')          
          NEUMANN(P, X, Y, Z, h, nx, ny, nz, sit);             
}


void DIRICHLET(double ***P, double X, double Y, double Z, double h, int nx, int ny, int nz, int sit)
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
                    P[nx-1][j][k] = 2*(2*X*X-5*y*y+4*z*z-1/3) - P[nx-2][j][k];
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
                    P[i][ny-1][k] = 2*(2*x*x+4*z*z-5*Y*Y-1/3) - P[i][ny-2][k];
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
                    P[i][j][nz-1] = 2*(2*x*x-5*y*y+4*Z*Z-1/3) - P[i][j][nz-2];
               }
               else 
               {
              		P[i][j][0] = -P[i][j][1];
               		P[i][j][nz-1] = -P[i][j][nz-2];
               }
	      }
     } 
}


void NEUMANN(double ***P, double X, double Y, double Z, double h, int nx, int ny, int nz, int sit)
{        
     int i, j, k;
     
     double x, y, z;
          
     if(sit == 1)          
          CELL_REFERENCE(P, nx, ny, nz);
     
     for(j=0; j<ny; j++)
     {
          y = 0.5*h + (j-1)*h;  
                      
	      for(k=0; k<nz; k++)
	      {
	           z = 0.5*h + (k-1)*h;
	           
               if(sit == 1)     
               {               
                    P[0][j][k] = P[1][j][k]; 
                    P[nx-1][j][k] = P[nx-2][j][k] + h*(4*X);
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
                    P[i][ny-1][k] = P[i][ny-2][k] + h*(-10*Y);
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
                    P[i][j][nz-1] = P[i][j][nz-2] + h*(8*Z);
               }
               else 
               {
              		P[i][j][0] = P[i][j][1];
               		P[i][j][nz-1] = P[i][j][nz-2];
               }
	      }
     } 
}


