//SISTEMA

#include <cstdlib>
#include <iostream>
using std::cout;
using std::cin;
using std::endl;
#include <math.h>
#include <stdio.h>

#include "auxiliares.h"
#include "propriedades.h"
#include "sistema.h"
#include "solvers.h"
#include "utilitarios.h"


/*

Boundary Condition

W : Dirichlet (Tw)
E : Dirichlet (Te)
S : Homogeneous Neumann
N : Homogeneous Neumann

*/

void coefficients(char **flag, double ***Ap, double ***Ae, double ***Aw, double ***An, double ***As, double **B, double ***Cp, double ***K, double **T0, double **Cp0, double **f, double **f0, int *nx, int *ny, double *h, double ro, double L, double dt, int levels, double Tw, double Te) 
{
     
     int i, j, l;
     
     double Kw, Ke, Kn, Ks;
     
     for(l=levels-1; l>=0; l--)
     {
          if(l == (levels-1))
          {               
               for(i=1; i<nx[l]-1; i++)
                    for(j=1; j<ny[l]-1; j++)
                    {                             
                         Kw = 0.5*(K[l][i-1][j]+K[l][i][j]);
                         Ke = 0.5*(K[l][i+1][j]+K[l][i][j]);
                         Kn = 0.5*(K[l][i][j+1]+K[l][i][j]);   
                         Ks = 0.5*(K[l][i][j-1]+K[l][i][j]);
                             
                         if(flag[i][j] == 'P')
                         {
                              Ae[l][i][j] = -Ke/(h[l]*h[l]);
                              Aw[l][i][j] = -Kw/(h[l]*h[l]);
                              An[l][i][j] = -Kn/(h[l]*h[l]);
                              As[l][i][j] = -Ks/(h[l]*h[l]);
                              Ap[l][i][j] = (Ke+Kw+Kn+Ks)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]);
                         }
                         else if(flag[i][j] == 'W')
                         {
                              Ae[l][i][j] = -Ke/(h[l]*h[l]);
                              Aw[l][i][j] = 0;
                              An[l][i][j] = -Kn/(h[l]*h[l]);
                              As[l][i][j] = -Ks/(h[l]*h[l]);
                              Ap[l][i][j] = (Ke+2*Kw+Kn+Ks)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]) + 2*Tw*Kw/(h[l]*h[l]);
                         }
                         else if(flag[i][j] == 'E')
                         {
                              Ae[l][i][j] = 0;
                              Aw[l][i][j] = -Kw/(h[l]*h[l]);
                              An[l][i][j] = -Kn/(h[l]*h[l]);
                              As[l][i][j] = -Ks/(h[l]*h[l]);
                              Ap[l][i][j] = (2*Ke+Kw+Kn+Ks)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]) + 2*Te*Ke/(h[l]*h[l]);
                         }                         
                         else if(flag[i][j] == 'N')
                         {
                              Ae[l][i][j] = -Ke/(h[l]*h[l]);
                              Aw[l][i][j] = -Kw/(h[l]*h[l]);
                              An[l][i][j] = 0;
                              As[l][i][j] = -Ks/(h[l]*h[l]);
                              Ap[l][i][j] = (Ke+Kw+Ks)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]);
                         }
                         else if(flag[i][j] == 'S')
                         {
                              Ae[l][i][j] = -Ke/(h[l]*h[l]);
                              Aw[l][i][j] = -Kw/(h[l]*h[l]);
                              An[l][i][j] = -Kn/(h[l]*h[l]);
                              As[l][i][j] = 0;
                              Ap[l][i][j] = (Ke+Kw+Kn)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]);
                         }
                         else if(flag[i][j] == 'A')
                         {
                              Ae[l][i][j] = -Ke/(h[l]*h[l]);
                              Aw[l][i][j] = 0;
                              An[l][i][j] = 0;
                              As[l][i][j] = -Ks/(h[l]*h[l]);
                              Ap[l][i][j] = (Ke+2*Kw+Ks)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]) + 2*Tw*Kw/(h[l]*h[l]);
                         }
                         else if(flag[i][j] == 'B')
                         {
                              Ae[l][i][j] = 0;
                              Aw[l][i][j] = -Kw/(h[l]*h[l]);
                              An[l][i][j] = 0;
                              As[l][i][j] = -Ks/(h[l]*h[l]);
                              Ap[l][i][j] = (2*Ke+Kw+Ks)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]) + 2*Te*Ke/(h[l]*h[l]);
                         }
                         else if(flag[i][j] == 'C')
                         {
                              Ae[l][i][j] = -Ke/(h[l]*h[l]);
                              Aw[l][i][j] = 0;
                              An[l][i][j] = -Kn/(h[l]*h[l]);
                              As[l][i][j] = 0;
                              Ap[l][i][j] = (Ke+2*Kw+Kn)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]) + 2*Tw*Kw/(h[l]*h[l]);
                         }
                         else if(flag[i][j] == 'D')
                         {
                              Ae[l][i][j] = 0;
                              Aw[l][i][j] = -Kw/(h[l]*h[l]);
                              An[l][i][j] = -Kn/(h[l]*h[l]);
                              As[l][i][j] = 0;
                              Ap[l][i][j] = (2*Ke+Kw+Kn)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]) + 2*Te*Ke/(h[l]*h[l]);
                         }                             
                    }                    
          }
          
          else if(l != (levels-1))
          {
               for(i=1; i<nx[l]-1; i++)
                    for(j=1; j<ny[l]-1; j++)
                    {                             
                         Kw = 0.5*(K[l][i-1][j]+K[l][i][j]);
                         Ke = 0.5*(K[l][i+1][j]+K[l][i][j]);
                         Kn = 0.5*(K[l][i][j+1]+K[l][i][j]);   
                         Ks = 0.5*(K[l][i][j-1]+K[l][i][j]);
                             
                         if(flag[i][j] == 'P')
                         {
                              Ae[l][i][j] = -Ke/(h[l]*h[l]);
                              Aw[l][i][j] = -Kw/(h[l]*h[l]);
                              An[l][i][j] = -Kn/(h[l]*h[l]);
                              As[l][i][j] = -Ks/(h[l]*h[l]);
                              Ap[l][i][j] = (Ke+Kw+Kn+Ks)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]);
                         }
                         else if(flag[i][j] == 'W')
                         {
                              Ae[l][i][j] = -Ke/(h[l]*h[l]);
                              Aw[l][i][j] = 0;
                              An[l][i][j] = -Kn/(h[l]*h[l]);
                              As[l][i][j] = -Ks/(h[l]*h[l]);
                              Ap[l][i][j] = (Ke+2*Kw+Kn+Ks)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]);
                         }
                         else if(flag[i][j] == 'E')
                         {
                              Ae[l][i][j] = 0;
                              Aw[l][i][j] = -Kw/(h[l]*h[l]);
                              An[l][i][j] = -Kn/(h[l]*h[l]);
                              As[l][i][j] = -Ks/(h[l]*h[l]);
                              Ap[l][i][j] = (2*Ke+Kw+Kn+Ks)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]);
                         }                         
                         else if(flag[i][j] == 'N')
                         {
                              Ae[l][i][j] = -Ke/(h[l]*h[l]);
                              Aw[l][i][j] = -Kw/(h[l]*h[l]);
                              An[l][i][j] = 0;
                              As[l][i][j] = -Ks/(h[l]*h[l]);
                              Ap[l][i][j] = (Ke+Kw+Ks)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]);
                         }
                         else if(flag[i][j] == 'S')
                         {
                              Ae[l][i][j] = -Ke/(h[l]*h[l]);
                              Aw[l][i][j] = -Kw/(h[l]*h[l]);
                              An[l][i][j] = -Kn/(h[l]*h[l]);
                              As[l][i][j] = 0;
                              Ap[l][i][j] = (Ke+Kw+Kn)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]);
                         }
                         else if(flag[i][j] == 'A')
                         {
                              Ae[l][i][j] = -Ke/(h[l]*h[l]);
                              Aw[l][i][j] = 0;
                              An[l][i][j] = 0;
                              As[l][i][j] = -Ks/(h[l]*h[l]);
                              Ap[l][i][j] = (Ke+2*Kw+Ks)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]);
                         }
                         else if(flag[i][j] == 'B')
                         {
                              Ae[l][i][j] = 0;
                              Aw[l][i][j] = -Kw/(h[l]*h[l]);
                              An[l][i][j] = 0;
                              As[l][i][j] = -Ks/(h[l]*h[l]);
                              Ap[l][i][j] = (2*Ke+Kw+Ks)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]);
                         }
                         else if(flag[i][j] == 'C')
                         {
                              Ae[l][i][j] = -Ke/(h[l]*h[l]);
                              Aw[l][i][j] = 0;
                              An[l][i][j] = -Kn/(h[l]*h[l]);
                              As[l][i][j] = 0;
                              Ap[l][i][j] = (Ke+2*Kw+Kn)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]);
                         }
                         else if(flag[i][j] == 'D')
                         {
                              Ae[l][i][j] = 0;
                              Aw[l][i][j] = -Kw/(h[l]*h[l]);
                              An[l][i][j] = -Kn/(h[l]*h[l]);
                              As[l][i][j] = 0;
                              Ap[l][i][j] = (2*Ke+Kw+Kn)/(h[l]*h[l]) + (ro/dt)*Cp[l][i][j];
                              B[i][j] = (ro/dt)*Cp0[i][j]*T0[i][j] - (ro*L/dt)*(f[i][j]-f0[i][j]);
                         }                             
                    } 
          } 
     }    
}              


void residuo(double **R, double **T, double **rhs, double **Ap, double **Ae, double **Aw, double **An, double **As, int nx, int ny)
{

        int i, j;

        for(i=1; i<(nx-1); i++)
                for(j=1; j<(ny-1); j++)
                        R[i][j] = rhs[i][j] - (Ae[i][j]*T[i+1][j] + Aw[i][j]*T[i-1][j] + An[i][j]*T[i][j+1] + As[i][j]*T[i][j-1] + Ap[i][j]*T[i][j]);

}
     
