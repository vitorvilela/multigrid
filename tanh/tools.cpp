//tools.cpp

#include <cstdlib>
#include <iostream>
using std::cout;

#include <stdio.h>
#include <math.h>

#include "system.h"
#include "solvers.h"
#include "tools.h"
#include "exact.h"


double INFINITY_NORM(double ***Pex, double ***P, double ***E, int nx, int ny, int nz)
{
     int i, j, k;
     
     double MAX = 0.0;
        
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               for(k=1; k<(nz-1); k++)
               {
         	        E[i][j][k] = Pex[i][j][k] - P[i][j][k];
         	        
 	                if(fabs(E[i][j][k]) > MAX)
                         MAX = fabs(E[i][j][k]);
               }
               
     return MAX;     
}


double MAXIMUM(double ***R, int nx, int ny, int nz)
{
     int i, j, k;

     double MAX = 0.0;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
		       for(k=1; k<(nz-1); k++)
               {
                    if(fabs(R[i][j][k]) > MAX)
                         MAX = fabs(R[i][j][k]);
               }
               
     return MAX;
}


double L2_NORM(double ***E, double h, int nx, int ny, int nz)
{
     int i, j, k;
     
	 double L2 = 0.0;

	 for(i=1; i<(nx-1); i++)
	      for(j=1; j<(ny-1); j++)
	           for(k=1; k<(nz-1); k++)
			        L2 += pow(E[i][j][k],2.0);    

     L2 = sqrt(L2*pow(h,3.0));

	 return L2;
}


void ZERO(double ****e, int *nx, int *ny, int *nz, int levels)
{
     int t, i, j, k;

     for(t=0; t<levels; t++)
          for(i=0; i<nx[t]; i++)
               for(j=0; j<ny[t]; j++)
                    for(k=0; k<nz[t]; k++)
                         e[t][i][j][k] = 0;
}


void CORRECT(double ***P, double ***e, int nx, int ny, int nz)
{
     int i, j, k;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
		       for(k=1; k<(nz-1); k++)
                    P[i][j][k] += e[i][j][k];
}


void INTEGRATION_REFERENCE(double ***P, double h, int nx, int ny, int nz)
{
     int i, j, k;
     
	 double SUM = 0.0;	 

	 for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               for(k=1; k<(nz-1); k++)
			        SUM = SUM + P[i][j][k];
			        
     SUM = SUM*h*h*h;
	 
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               for(k=1; k<(nz-1); k++)
                    P[i][j][k] = P[i][j][k] - SUM;	    
}


void CELL_REFERENCE(double ***P, int nx, int ny, int nz)
{
     int i, j, k;
     
     double CELL = P[1][1][1];
     
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               for(k=1; k<(nz-1); k++)
			        P[i][j][k] = P[i][j][k] - CELL;
}


void EQUAL(double ***rhs, double ***R, int nx, int ny, int nz)
{
     int i, j, k;

     for(i=0; i<nx; i++)
          for(j=0; j<ny; j++)
               for(k=0; k<nz; k++)
                    rhs[i][j][k] = R[i][j][k];
}

double INTERNAL_PRODUCT(double ***M1, double ***M2, int nx, int ny, int nz)
{
     int i, j, k;
     
	 double SUM = 0.0;

	 for(i=1; i<(nx-1); i++)
	     for(j=1; j<(ny-1); j++)
	          for(k=1; k<(nz-1); k++)
			       SUM += M1[i][j][k]*M2[i][j][k];

	 return SUM;
}


void INTERPOLATION(double ****e, int *nx, int *ny, int *nz, int t)
{     
     int i, j, k;

     int limitX = nx[t]-2;
	 int limitY = ny[t]-2;
	 int limitZ = nz[t]-2;

     for(i=1; i<=limitX; i++)
          for(j=1; j<=limitY; j++)
               for(k=1; k<=limitZ; k++)
               {                    
                    e[t+1][2*i-1][2*j-1][2*k-1] = (1.0/64.0)*e[t][i-1][j-1][k-1] + (3.0/64.0)*e[t][i][j-1][k-1] + (9.0/64.0)*e[t][i][j][k-1] + (3.0/64.0)*e[t][i-1][j][k-1] + (3.0/64.0)*e[t][i-1][j-1][k] + (9.0/64.0)*e[t][i][j-1][k] + (27.0/64.0)*e[t][i][j][k] + (9.0/64.0)*e[t][i-1][j][k];    
                    
                    e[t+1][2*i-1][2*j-1][2*k] = (3.0/64.0)*e[t][i-1][j-1][k] + (9.0/64.0)*e[t][i][j-1][k] + (27.0/64.0)*e[t][i][j][k] + (9.0/64.0)*e[t][i-1][j][k] + (1.0/64.0)*e[t][i-1][j-1][k+1] + (3.0/64.0)*e[t][i][j-1][k+1] + (9.0/64.0)*e[t][i][j][k+1] + (3.0/64.0)*e[t][i-1][j][k+1];
                        
                    e[t+1][2*i-1][2*j][2*k-1] = (3.0/64.0)*e[t][i-1][j][k-1] + (9.0/64.0)*e[t][i][j][k-1] + (3.0/64.0)*e[t][i][j+1][k-1] + (1.0/64.0)*e[t][i-1][j+1][k-1] + (9.0/64.0)*e[t][i-1][j][k] + (27.0/64.0)*e[t][i][j][k] + (9.0/64.0)*e[t][i][j+1][k] + (3.0/64.0)*e[t][i-1][j+1][k];
                    
                    e[t+1][2*i-1][2*j][2*k] = (9.0/64.0)*e[t][i-1][j][k] + (27.0/64.0)*e[t][i][j][k] + (9.0/64.0)*e[t][i][j+1][k] + (3.0/64.0)*e[t][i-1][j+1][k] + (3.0/64.0)*e[t][i-1][j][k+1] + (9.0/64.0)*e[t][i][j][k+1] + (3.0/64.0)*e[t][i][j+1][k+1] + (1.0/64.0)*e[t][i-1][j+1][k+1];
		      
		            e[t+1][2*i][2*j-1][2*k-1] = (3.0/64.0)*e[t][i][j-1][k-1] + (1.0/64.0)*e[t][i+1][j-1][k-1] + (3.0/64.0)*e[t][i+1][j][k-1] + (9.0/64.0)*e[t][i][j][k-1] + (9.0/64.0)*e[t][i][j-1][k] + (3.0/64.0)*e[t][i+1][j-1][k] + (9.0/64.0)*e[t][i+1][j][k] + (27.0/64.0)*e[t][i][j][k];
		            
		            e[t+1][2*i][2*j-1][2*k] = (9.0/64.0)*e[t][i][j-1][k] + (3.0/64.0)*e[t][i+1][j-1][k] + (9.0/64.0)*e[t][i+1][j][k] + (27.0/64.0)*e[t][i][j][k] + (3.0/64.0)*e[t][i][j-1][k+1] + (1.0/64.0)*e[t][i+1][j-1][k+1] + (3.0/64.0)*e[t][i+1][j][k+1] + (9.0/64.0)*e[t][i][j][k+1];
		            
		            e[t+1][2*i][2*j][2*k-1] = (9.0/64.0)*e[t][i][j][k-1] + (3.0/64.0)*e[t][i+1][j][k-1] + (1.0/64.0)*e[t][i+1][j+1][k-1] + (3.0/64.0)*e[t][i][j+1][k-1] + (27.0/64.0)*e[t][i][j][k] + (9.0/64.0)*e[t][i+1][j][k] + (3.0/64.0)*e[t][i+1][j+1][k] + (9.0/64.0)*e[t][i][j+1][k];
		            
		            e[t+1][2*i][2*j][2*k] = (27.0/64.0)*e[t][i][j][k] + (9.0/64.0)*e[t][i+1][j][k] + (3.0/64.0)*e[t][i+1][j+1][k] + (9.0/64.0)*e[t][i][j+1][k] + (9.0/64.0)*e[t][i][j][k+1] + (3.0/64.0)*e[t][i+1][j][k+1] + (1.0/64.0)*e[t][i+1][j+1][k+1] + (3.0/64.0)*e[t][i][j+1][k+1];		            		            
               } 
}


void RESTRICTION(double ***R, double ****rhs, int *nx, int *ny, int *nz, int t)
{
	 int i, j, k;
	
	 int limitX = nx[t-1]-2;
     int limitY = ny[t-1]-2;
     int limitZ = nz[t-1]-2;

     for(i=1; i<=limitX; i++)
          for(j=1; j<=limitY; j++)
               for(k=1; k<=limitZ; k++)
                    rhs[t-1][i][j][k] = 0.125*(R[2*i-1][2*j-1][2*k-1] + R[2*i-1][2*j-1][2*k] + R[2*i-1][2*j][2*k-1] + R[2*i-1][2*j][2*k] + R[2*i][2*j-1][2*k-1] + R[2*i][2*j-1][2*k] + R[2*i][2*j][2*k-1] + R[2*i][2*j][2*k]);
	
}




