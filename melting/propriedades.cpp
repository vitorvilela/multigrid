//PROPRIEDADES

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



void set_K(double ***K, double **T, int *nx, int *ny, int niveis)
{
     int t, i, j; 
     
     int pausa;
     
     t = niveis-1;                   
                     
     double Ts = 1400;   
        
     for(i=1; i<(nx[t]-1); i++)
          for(j=1; j<(ny[t]-1); j++)
          {
               if(T[i][j] <= Ts)
                    K[t][i][j] = 14.3 + 0.01983*T[i][j] - 5.451E-06*T[i][j]*T[i][j];
               else
                    K[t][i][j] = 31.37804;
          } 
     
     Homogeneous_Neumann(K[t], nx[t], ny[t]);
               
     for(t=(niveis-1); t>0; t--)
     {
          restringe(K[t], K, nx, ny, t);
          Homogeneous_Neumann(K[t-1], nx[t-1], ny[t-1]);
     }
     
}

void set_Cp0(double **Cp0, double **T, int nx, int ny)
{
     int i, j;
          
     double Ts = 1400;    
     
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
          {
               if(T[i][j] <= Ts)
                    Cp0[i][j] = 460.5 + 0.4257*T[i][j] - 5.05E-04*T[i][j]*T[i][j] + 2.6608E-07*T[i][j]*T[i][j]*T[i][j];
               else
                    Cp0[i][j] = 796.584;
          }
}

void set_Cp(double ***Cp, double **Cp0, double **T, double **Told, int *nx, int *ny, int niveis)
{
     int t, i, j;
     
     t = niveis-1;
     
     double Ts = 1400;
          
     for(i=1; i<(nx[t]-1); i++)
          for(j=1; j<(ny[t]-1); j++)
          {
               if(T[i][j] <= Ts)
                    Cp[t][i][j] = Cp0[i][j] + (0.4257 - 1.01E-03*T[i][j] + 7.9824E-07*T[i][j]*T[i][j])*(T[i][j]-Told[i][j]);
               else
                    Cp[t][i][j] = Cp0[i][j];              
          }
          
     for(t=(niveis-1); t>0; t--)     
          restringe(Cp[t], Cp, nx, ny, t);         
}
    
void set_f0(double **f0, double **T, int nx, int ny)
{
     int i, j;
     
     double Ts = 1400;
     double Tl = 1455;    
     
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
          {
               if(T[i][j] < Ts)
                    f0[i][j] = 0;          
               else if((T[i][j] >= Ts) && (T[i][j] <= Tl))
                    f0[i][j] = (T[i][j]-Ts)/(Tl-Ts);
               else
                    f0[i][j] = 1;
          }   
}

void set_f(double **f, double **f0, double **T, double **Told, int nx, int ny)
{
     int i, j;
     
     double Ts = 1400;
     double Tl = 1455; 
     
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
          {
               if((T[i][j] >= Ts) && (T[i][j] <= Tl))
               {
                     f[i][j] = f0[i][j] + (1/(Tl-Ts))*(T[i][j]-Told[i][j]);
                    
                     if(f[i][j] > 1)   f[i][j] = 1.0;
                     else if(f[i][j] < 0)   f[i][j] = 0.0;   
               }                    
               else
                    f[i][j] = f0[i][j];
          }     
}    
 
   
