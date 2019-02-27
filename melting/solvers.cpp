//SOLVERS

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



void redblack(double **T, double **B, double **Ap, double **Ae, double **Aw, double **An, double **As, int nx, int ny, double L, double D, double h, double w, int itr, int sit, int fun, char fron)
{
     
     int i, j, irb, jrb, k;
     int it;
     
     

	 for(it=1; it<=itr; it++)
     {
	      for(k=-1; k<=1; k+=2)
   		  {
               jrb = k;
     		   for(j=1; j<(ny-1); j++)
     		   {
                    jrb = -jrb;
                    irb = (jrb+1)/2;
     		        for(i=(irb+1); i<(nx-1); i+=2)
             		     T[i][j] = (w/Ap[i][j])*(B[i][j]-(Ae[i][j]*T[i+1][j]+Aw[i][j]*T[i-1][j]+An[i][j]*T[i][j+1]+As[i][j]*T[i][j-1]))+(1-w)*T[i][j];
        		}
                
                
		  }
	 }
}


void sor(double **T, double **K, double **B, double **Ap, double **Ae, double **Aw, double **An, double **As, int nx, int ny, double L, double D, double h, double w, int itr)
{
     int i, j;
     int it;

     for(it=1; it<=itr; it++)                   
          for(i=1; i<(nx-1); i++)
               for(j=1; j<(ny-1); j++)
                    T[i][j] = (w/Ap[i][j])*(B[i][j]-(Ae[i][j]*T[i+1][j]+Aw[i][j]*T[i-1][j]+An[i][j]*T[i][j+1]+As[i][j]*T[i][j-1]))+(1-w)*T[i][j];
}

/*
void jacobi(double **P, double **B, double **Ap, int nx, int ny)
{
     int i, j;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               P[i][j] = (B[i][j]/Ap[i][j]);
}

void cg(double **R1, double **R2, double **z1, double **z2, double **p, double **W, double **P, double **B, double ***Ap, double ***Ae, double ***Aw, double ***As, double ***An, int *nx, int *ny, double L, double D, double *h, int fun, char fron, int niveis, double eps, char pre, int *itr_u, int *itr_d, double *w, double ***e, double ***rhs, double ***R, double ***Zero)
{
     cout << "\n\n            CG_PRECONDICIONADO            \n\n";  
     
     int i, j;
     int t = niveis-1;
     
     int cont = 0;

     double alfa, beta;
     double rsd;
     
     int pausa;     
  
  
     fronteira(P, h[t], nx[t], ny[t], L, D, 1, fun, fron);
     residuo(R1, P, B, Ap[t], Ae[t], Aw[t], An[t], As[t], nx[t], ny[t]);
     rsd = maximo(R1, nx[t], ny[t]);
     cout << "\nRESIDUO INICIAL = " << rsd << endl;

     if(pre == 's')
          igual(z1, R1, nx[t], ny[t]);
     else if(pre == 'j')     
          jacobi(z1, R1, Ap[t], nx[t], ny[t]); 
     else if(pre == 'm')
          multigrid(e, R, rhs, z1, R1, Ap, Ae, Aw, Zero, Zero, itr_u, itr_d, nx, ny, L, D, h, w, fun, fron, eps, niveis);
     
              
     igual(p, z1, nx[t], ny[t]);  
     
     if(fron == 'd' || fron == 'p')
          fronteira(p, h[t], nx[t], ny[t], L, D, 0, fun, fron);
     else
          if(t==(niveis-1)) 
               fronteira(p, h[t], nx[t], ny[t], L, D, 0, fun, fron);      
               
     for(i=1; i<(nx[t]-1); i++)
         for(j=1; j<(ny[t]-1); j++)
              W[i][j] = Ap[t][i][j]*p[i][j] + Ae[t][i][j]*p[i+1][j] + Aw[t][i][j]*p[i-1][j] + An[t][i][j]*p[i][j+1] + As[t][i][j]*p[i][j-1];
     
     alfa = escalar(R1, z1, nx[t], ny[t])/escalar(p, W, nx[t], ny[t]);
          
     do
     {   
          cont++;         
         
          for(i=1; i<(nx[t]-1); i++)
               for(j=1; j<(ny[t]-1); j++)
                    P[i][j] += alfa*p[i][j];          
         
          fronteira(P, h[t], nx[t], ny[t], L, D, 1, fun, fron);     

          residuo(R2, P, B, Ap[t], Ae[t], Aw[t], An[t], As[t], nx[t], ny[t]);
          rsd = maximo(R2, nx[t], ny[t]);
          cout << "CONT: " << cont << "      RESIDUO = " << rsd << endl;
          
          if(pre == 's')
               igual(z2, R2, nx[t], ny[t]);
          else if(pre == 'j')     
               jacobi(z2, R2, Ap[t], nx[t], ny[t]); 
          else if(pre == 'm')
               multigrid(e, R, rhs, z2, R2, Ap, Ae, Aw, Zero, Zero, itr_u, itr_d, nx, ny, L, D, h, w, fun, fron, eps, niveis);
                          
                       
          beta = escalar(R2, z2, nx[t], ny[t])/escalar(R1, z1, nx[t], ny[t]);
          
          for(i=1; i<(nx[t]-1); i++)
               for(j=1; j<(ny[t]-1); j++)
                    p[i][j] = z2[i][j] + beta*p[i][j];          
          
          if(fron == 'd' || fron == 'p')
               fronteira(p, h[t], nx[t], ny[t], L, D, 0, fun, fron);
          else
               if(t==(niveis-1)) 
                    fronteira(p, h[t], nx[t], ny[t], L, D, 0, fun, fron);
          
          for(i=1; i<(nx[t]-1); i++)
               for(j=1; j<(ny[t]-1); j++)
                    W[i][j] = Ap[t][i][j]*p[i][j] + Ae[t][i][j]*p[i+1][j] + Aw[t][i][j]*p[i-1][j] + An[t][i][j]*p[i][j+1] + As[t][i][j]*p[i][j-1];
                    
          alfa = escalar(R2, z2, nx[t], ny[t])/escalar(p, W, nx[t], ny[t]);
          
          igual(R1, R2, nx[t], ny[t]);
          
          igual(z1, z2, nx[t], ny[t]);
          
          //cin >> pausa;                     
     }
     while(rsd > eps);   
                      
}
*/
/*
void multigrid(double ***e, double ***R, double ***rhs, double **T, double **B, double ***Ap, double ***Ae, double ***Aw, double ***An, double ***As, int *itr_u, int *itr_d, int *nx, int *ny, double L, double D, double *h, double *w, int fun, char fron, double eps, int niveis)
{
     //cout << "\n\n         MULTIGRID\n\n";
     
     int pausa;
     
     int t, i, j;
     int NC = 0;
     double rsd = 1;
    
     
     do
     {
          t = niveis-1;
          
          NC++;          
                
          redblack(T, B, Ap[t], Ae[t], Aw[t], An[t], As[t], nx[t], ny[t], L, D, h[t], w[t], itr_d[t], 1, fun, fron);           
          residuo(R[t], T, B, Ap[t], Ae[t], Aw[t], An[t], As[t], nx[t], ny[t]);		                
          rsd = maximo(R[t], nx[t], ny[t]);
          cout << "NC-down: " << NC << " --> rsd: " << rsd << endl;

          while(t > 0)
          {               
               restringe(R[t], rhs, nx, ny, t);
               t--;                       
               redblack(e[t], rhs[t], Ap[t], Ae[t], Aw[t], An[t], As[t], nx[t], ny[t], L, D, h[t], w[t], itr_d[t], 0, fun, fron);               	
               residuo(R[t], e[t], rhs[t], Ap[t], Ae[t], Aw[t], An[t], As[t], nx[t], ny[t]);		  			   
           }

           while(t < (niveis-1))
           {                                 
                interpola(e, nx, ny, t);                        
                t++;
                		
                if(t != (niveis-1))							
                     redblack(e[t], rhs[t], Ap[t], Ae[t], Aw[t], An[t], As[t], nx[t], ny[t], L, D, h[t], w[t], itr_u[t], 0, fun, fron);
   			                                    
                else if(t == (niveis-1))
                {
                     corrigeMAIS(P, e[t], nx[t], ny[t]);
                     zero(e, nx, ny, niveis);                          
                }
           }
           
           //cin >> pausa;
           
     }
     while(rsd > eps);

	 cout << "\n";
}
*/
/*

void nested(double ***e, double ***R, double ***rhs, double ***Ap, double ***Ae, double ***Aw, double ***As, double ***An, double ***Pex, double ***B, double ***P, int *nx, int *ny, double *h, int niveis, int fun, char fron, double L, double D, int *itr_u, int *itr_d, double eps, double *w)
{
      cout << "\n\n                  NESTED\n\n";
     
      int t, i, j;
            
      int pausa;
      
      int NIVEIS;
         
      for(t=0; t<niveis; t++)
      {
           cout << "\n\nT = " << t << "\n\n";                
         
           NIVEIS = t + 1;              
           multigrid(e, R, rhs, P[t], B[t], Ap, Ae, Aw, An, As, itr_u, itr_d, nx, ny, L, D, h, w, fun, fron, eps, NIVEIS);
                          
          
           if(t != (niveis-1))
                interpola(P, nx, ny, t); 
                
           //cin >> pausa;
      }   
}




void MGgM(double ***p, double ***R, double ***K, double ***S, double ***Ap, double ***Ae, double ***Aw, double ***As, double ***An, double ***B, double ***P, int *nx, int *ny, double *h, int niveis, int fun, char fron, double L, double D, double eps, char pre)
{     
     cout << "\n\n                          Multigrid Gradient Algorithm\n\n";     
     
     int t, i, j;     
     int pausa; 
     int NC = 0;
     
     double rsd;
     double epR;
     double seKp; 
     double eKp;    
          
     t = niveis-1;     
     fronteira(P[t], h[t], nx[t], ny[t], L, D, 1, fun, fron);     
     residuo(R[t], P[t], B[t], Ap[t], Ae[t], Aw[t], An[t], As[t], nx[t], ny[t]);     
     rsd = maximo(R[t], nx[t], ny[t]);
     cout << "NC = " << NC << " - RESIDUO INICIAL = " << rsd << endl; 
     
     if(pre == 's')     
          igual(p[t], R[t], nx[t], ny[t]);
     else if(pre == 'j')
          jacobi(p[t], R[t], Ap[t], nx[t], ny[t]);         
         
     do
     {      
          NC++;     
               
          while(t>0)
          {
               restringe(R[t], R, nx, ny, t);
               t--;    
                          
               if(pre == 's')     
                    igual(p[t], R[t], nx[t], ny[t]);
               else if(pre == 'j')
                    jacobi(p[t], R[t], Ap[t], nx[t], ny[t]);
          }    
     
          t = 0;         
     
          while(t<=(niveis-1))
          {    
               if(fron == 'd' || fron == 'p')
                    fronteira(p[t], h[t], nx[t], ny[t], L, D, 0, fun, fron); 
               else
                   if(t==(niveis-1))                                             
                        fronteira(p[t], h[t], nx[t], ny[t], L, D, 0, fun, fron);                       
               
               for(i=1; i<(nx[t]-1); i++)
                    for(j=1; j<(ny[t]-1); j++)
                         K[t][i][j] = Ap[t][i][j]*p[t][i][j] + Ae[t][i][j]*p[t][i+1][j] + Aw[t][i][j]*p[t][i-1][j] + An[t][i][j]*p[t][i][j+1] + As[t][i][j]*p[t][i][j-1];
                                           
               for(int tt=t; tt>0; tt--)
                    restringe(K[tt], K, nx, ny, tt);              
                
               for(int tt=0; tt<niveis; tt++)
                    for(i=0; i<nx[tt]; i++)
                         for(j=0; j<ny[tt]; j++)
                             S[tt][i][j] = 0;               
                     
               for(int tt=0; tt<t; tt++)
               {                     
                    eKp = escalar(K[tt], p[tt], nx[tt], ny[tt]);
                                    
                    for(i=1; i<(nx[tt]-1); i++)
                         for(j=1; j<(ny[tt]-1); j++)   
                              S[tt][i][j] += p[tt][i][j]*eKp;                    
                             
                    fronteira(S[tt], h[tt], nx[tt], ny[tt], L, D, 0, fun, fron);                        
                    interpola(S, nx, ny, tt);                    
               }               
                          
               for(i=1; i<(nx[t]-1); i++)
                    for(j=1; j<(ny[t]-1); j++)
                         p[t][i][j] -= S[t][i][j];
                         
               if(fron == 'd' || fron == 'p')
                    fronteira(p[t], h[t], nx[t], ny[t], L, D, 0, fun, fron); 
               else
                   if(t==(niveis-1))                                             
                        fronteira(p[t], h[t], nx[t], ny[t], L, D, 0, fun, fron);
                                       
               for(i=1; i<(nx[t]-1); i++)
                    for(j=1; j<(ny[t]-1); j++)
                         K[t][i][j] = Ap[t][i][j]*p[t][i][j] + Ae[t][i][j]*p[t][i+1][j] + Aw[t][i][j]*p[t][i-1][j] + An[t][i][j]*p[t][i][j+1] + As[t][i][j]*p[t][i][j-1];
                                          
               seKp = sqrt(escalar(K[t], p[t], nx[t], ny[t]));
                                  
               for(i=1; i<(nx[t]-1); i++)
                    for(j=1; j<(ny[t]-1); j++)
                         p[t][i][j] = p[t][i][j]/seKp;
                     
               t++;               
          }
          
          for(int tt=0; tt<niveis; tt++)
               for(i=0; i<nx[tt]; i++)
                    for(j=0; j<ny[tt]; j++)
                         S[tt][i][j] = 0;      
     
          t = 0;
          
          epR = escalar(p[t], R[t], nx[t], ny[t]);          
     
          for(i=1; i<(nx[0]-1); i++)
               for(j=1; j<(ny[0]-1); j++)
                    S[0][i][j] = p[0][i][j]*epR;
          
          while(t<(niveis-1))
          {                                       
               fronteira(S[t], h[t], nx[t], ny[t], L, D, 0, fun, fron);                             
               interpola(S, nx, ny, t);  
                         
               t++;             
               
               epR = escalar(p[t], R[t], nx[t], ny[t]);
                
               for(i=1; i<(nx[t]-1); i++)
                    for(j=1; j<(ny[t]-1); j++)   
                         S[t][i][j] += p[t][i][j]*epR;                     
          }    
     
          t = niveis-1;    
          corrigeMAIS(P[t], S[t], nx[t], ny[t]);
          fronteira(P[t], h[t], nx[t], ny[t], L, D, 1, fun, fron);
          residuo(R[t], P[t], B[t], Ap[t], Ae[t], Aw[t], An[t], As[t], nx[t], ny[t]);              
          rsd = maximo(R[t], nx[t], ny[t]);
          cout << "\nNC = " << NC << " - RESIDUO = " << rsd << endl; 
                
          if(pre == 's')     
               igual(p[t], R[t], nx[t], ny[t]);
          else if(pre == 'j')
               jacobi(p[t], R[t], Ap[t], nx[t], ny[t]);            
     }
     while(rsd > eps);
     
}
*/

