//SOLVERS

#include <cstdlib>
#include <iostream>

using std::cout;
using std::cin;
using std::endl;

#include <stdio.h>
#include <math.h>

#include "sistema.h"
#include "solvers.h"
#include "utilitarios.h"
#include "exatas.h"
#include "fronteiras.h"

/*
void nested(double ****Pex, double ****E, double ****e, double ****R, double ****rhs, double ****z1, double ****z2, double ****R1, double ****R2, double ****p, double ****W, double ****Ap, double ****Ae, double ****Aw, double ****As, double ****An, double ****Au, double ****Ad, double ****Mp, double ****Me, double ****Mw, double ****Ms, double ****Mn, double ****Mu, double ****Md, double ****B, double ****P, int *nx, int *ny, int *nz, double *h, int niveis, int fun, char fron, double L, double D, double S, int *itr_u, int *itr_d, double eps, double *w)
{
      cout << "\n\n                  NESTED\n\n";
     
      int t, i, j, k;
            
      int pausa;
      
      int NIVEIS;
      double EPS;
         
      for(t=0; t<niveis; t++)
      {
           cout << "\n\nT = " << t << "\n\n";                
         
           NIVEIS = t + 1;
           EPS = eps;              
           //MG(Pex[t], E[t], NIVEIS, e, R, rhs, P[t], B[t], Ap, Ae, Aw, An, As, Au, Ad, itr_u, itr_d, nx, ny, nz, D, L, S, h, w, EPS, fun, fron);
           cg(Pex[t], P[t], B[t], E[t], z1[t], z2[t], R1[t], R2[t], p[t], W[t], R, rhs, e, Ap, Ae, Aw, An, As, Au, Ad, Mp, Me, Mw, Mn, Ms, Mu, Md, Ap, nx, ny, nz, D, L, S, h, w, itr_d, itr_u, eps, fun, fron, 'j', NIVEIS);

          
           if(t != (niveis-1))
                interpola3d(P, nx, ny, nz, t); 
                
           //cin >> pausa;
      }   
}
*/

int MG(int niveis, double ****e, double ****R, double ****rhs, double ***P, double ***B, double ****Ap, double ****Ae, double ****Aw, double ****An, double ****As, double ****Au, double ****Ad, int *itr_u, int *itr_d, int *nx, int *ny, int *nz, double D, double L, double S, double *h, double *w, double eps, char fron)
{
     
      cout << "\n\n                         MULTIGRID\n\n";
            
      int t, i, j, k;     
      int pausa;        
      int NC = 0;
      
	  double rsd, emax;
	
      do
      {             
           t = niveis-1;
           NC++;		
		
		   sor(P, B, Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t], D, L, S, h[t], w[t], itr_d[t], 1, fron);
           residuo(R[t], P, B, Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t]);
           rsd = maximo(R[t], nx[t], ny[t], nz[t]);
           cout << "NC = " << NC << "     :     RESIDUO = " << rsd << endl;
           
           while(t > 0)
           {
                restringe3d(R[t], rhs, nx, ny, nz, t);
                t--;
                sor(e[t], rhs[t], Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t], D, L, S, h[t], w[t], itr_d[t], 0, fron);	
			    residuo(R[t], e[t], rhs[t], Ap[t], Ae[t], Aw[t], An[t], As[t],Au[t], Ad[t], nx[t], ny[t], nz[t]);
           }

           while(t < (niveis-1))
           {
                interpola3d(e, nx, ny, nz, t);
                t++;
			
			    if(t != (niveis-1))
				     sor(e[t], rhs[t], Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t], D, L, S, h[t], w[t], itr_u[t], 0, fron);
			         			
                else if(t == (niveis-1))
                {
                     corrige(P, e[t], nx[t], ny[t], nz[t]);
                     zero(e, nx, ny, nz, niveis);				     
                }
            }
            
            //cin >> pausa;
                                   
     }
     while(rsd > eps);
     
     return NC;
}


void sor(double ***Pe, double ***B, double ***Ap, double ***Ae, double ***Aw, double ***An, double ***As, double ***Au, double ***Ad, int nx, int ny, int nz, double D, double L, double S, double h, double w, int itr, int sit, char fron)
{
     int it, i, j, k;
     
     Fronteira_P(Pe, h, nx, ny, nz, fron, sit);
     //Cavidade_pressao(nx, ny, nz, Pe, h, sit);
     
     for(it=1; it<=itr; it++)
     {
  	      for(i=1; i<(nx-1); i++)
               for(j=1; j<(ny-1); j++)
			        for(k=1; k<(nz-1); k++)
                         Pe[i][j][k] = (w/Ap[i][j][k])*(B[i][j][k]-(Ae[i][j][k]*Pe[i+1][j][k]+
                                      Aw[i][j][k]*Pe[i-1][j][k]+An[i][j][k]*Pe[i][j+1][k]+
                                      As[i][j][k]*Pe[i][j-1][k]+Au[i][j][k]*Pe[i][j][k+1]+
                                      Ad[i][j][k]*Pe[i][j][k-1]))+(1-w)*Pe[i][j][k];
           
          Fronteira_P(Pe, h, nx, ny, nz, fron, sit);
          //Cavidade_pressao(nx, ny, nz, Pe, h, sit);
     }
}

/*
void jacobi(double ***P, double ***B, double ***Ap, int nx, int ny, int nz)
{
     int i, j, k;
               
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               for(k=1; k<(nz-1); k++)
                    P[i][j][k] = B[i][j][k]/Ap[i][j][k];    
}


void MGgM(double ***Pex, double ***E, double ****p, double ****R, double ****K, double ****S, double ****Ap, double ****Ae, double ****Aw, double ****As, double ****An, double ****Au, double ****Ad, double ****B, double ****P, int *nx, int *ny, int *nz, double *h, int niveis, int fun, char fron, double L, double D, double s, double eps, char pre)
{          
     cout << "\n\n                          Multigrid Gradient Algorithm\n\n"; 
     
     FILE *g1, *g2;
	 g1 = fopen("RSDxNC.dat", "w");	
     g2 = fopen("ERROxNC.dat", "w");
     fprintf(g1, "'NC'       'RSD'\n\n"); 
     fprintf(g2, "'NC'       'EMAX'\n\n");    
     
     int t, i, j, k;     
     int pausa; 
     int NC = 0;
     
     double rsd, emax;
     double epR;
     double seKp; 
     double eKp; 
     
     int it = 2;
     double w = 1;    
          
     t = niveis-1;    
     
     sor(P[t], B[t], Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t], D, L, s, h[t], w, it, 1, fun, fron);          
     fronteira(P[t], h[t], nx[t], ny[t], nz[t], D, L, s, fun, fron, 1);          
     residuo(R[t], P[t], B[t], Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t]);     
     rsd = maximo(R[t], nx[t], ny[t], nz[t]);
     cout << "NC = " << NC << " - RESIDUO INICIAL = " << rsd << endl; 
     
     fprintf(g1, "%i       %e\n", NC, rsd);     
     erro(Pex, P[t], E, nx[t], ny[t], nz[t]);
     emax = maximo(E, nx[t], ny[t], nz[t]);	      
     fprintf(g2, "%i       %e\n", NC, emax);
     
     
     if(pre == 's')     
          igual(p[t], R[t], nx[t], ny[t], nz[t]);
     else if(pre == 'j')
          jacobi(p[t], R[t], Ap[t], nx[t], ny[t], nz[t]);        
         
     do
     {           
          NC++;      
               
          while(t>0)
          {
               restringe3d(R[t], R, nx, ny, nz, t);
               t--;       
                       
               if(pre == 's')     
                    igual(p[t], R[t], nx[t], ny[t], nz[t]);
               else if(pre == 'j')
                    jacobi(p[t], R[t], Ap[t], nx[t], ny[t], nz[t]); 
          }         
     
          t = 0;       
     
          while(t<=(niveis-1))
          {    
               if(fron == 'd')                                            
                    fronteira(p[t], h[t], nx[t], ny[t], nz[t], D, L, s, fun, fron, 0);   
               else
                    if(t == (niveis-1))
                         fronteira(p[t], h[t], nx[t], ny[t], nz[t], D, L, s, fun, fron, 0);
                               
               for(i=1; i<(nx[t]-1); i++)
                    for(j=1; j<(ny[t]-1); j++)
                         for(k=1; k<(nz[t]-1); k++)
                              K[t][i][j][k] = Ap[t][i][j][k]*p[t][i][j][k] + Ae[t][i][j][k]*p[t][i+1][j][k] + Aw[t][i][j][k]*p[t][i-1][j][k] + An[t][i][j][k]*p[t][i][j+1][k] + As[t][i][j][k]*p[t][i][j-1][k] + Au[t][i][j][k]*p[t][i][j][k+1] + Ad[t][i][j][k]*p[t][i][j][k-1];
                                           
               for(int tt=t; tt>0; tt--)
                    restringe3d(K[tt], K, nx, ny, nz, tt);              
                
               for(i=1; i<(nx[0]-1); i++)
                    for(j=1; j<(ny[0]-1); j++)
                         for(k=1; k<(nz[0]-1); k++)
                              S[0][i][j][k] = 0;               
                     
               for(int tt=0; tt<t; tt++)
               {                  
                    eKp = escalar(K[tt], p[tt], nx[tt], ny[tt], nz[tt]);
                                    
                    for(i=1; i<(nx[tt]-1); i++)
                         for(j=1; j<(ny[tt]-1); j++)
                              for(k=1; k<(nz[tt]-1); k++)   
                                   S[tt][i][j][k] += p[tt][i][j][k]*eKp;
                          
                    for(i=1; i<(nx[tt+1]-1); i++)
                         for(j=1; j<(ny[tt+1]-1); j++)
                              for(k=1; k<(nz[tt+1]-1); k++)
                                   S[tt+1][i][j][k] = 0;
                             
                    fronteira(S[t], h[t], nx[t], ny[t], nz[t], D, L, s, fun, fron, 0);                        
                    interpola3d(S, nx, ny, nz, tt);                    
               }               
                          
               for(i=1; i<(nx[t]-1); i++)
                    for(j=1; j<(ny[t]-1); j++)
                         for(k=1; k<(nz[t]-1); k++)
                              p[t][i][j][k] -= S[t][i][j][k];
                         
               if(fron == 'd')                                            
                    fronteira(p[t], h[t], nx[t], ny[t], nz[t], D, L, s, fun, fron, 0);   
               else
                    if(t == (niveis-1))  
                         fronteira(p[t], h[t], nx[t], ny[t], nz[t], D, L, s, fun, fron, 0);
                                       
               for(i=1; i<(nx[t]-1); i++)
                    for(j=1; j<(ny[t]-1); j++)
                         for(k=1; k<(nz[t]-1); k++)
                              K[t][i][j][k] = Ap[t][i][j][k]*p[t][i][j][k] + Ae[t][i][j][k]*p[t][i+1][j][k] + Aw[t][i][j][k]*p[t][i-1][j][k] + An[t][i][j][k]*p[t][i][j+1][k] + As[t][i][j][k]*p[t][i][j-1][k] + Au[t][i][j][k]*p[t][i][j][k+1] + Ad[t][i][j][k]*p[t][i][j][k-1];
                                          
               seKp = sqrt(escalar(K[t], p[t], nx[t], ny[t], nz[t]));
                                  
               for(i=1; i<(nx[t]-1); i++)
                    for(j=1; j<(ny[t]-1); j++)
                         for(k=1; k<(nz[t]-1); k++)
                              p[t][i][j][k] = p[t][i][j][k]/seKp;
                     
               t++;               
          }    
               
          t = 0;    
          
          epR = escalar(p[t], R[t], nx[t], ny[t], nz[t]);          
     
          for(i=1; i<(nx[0]-1); i++)
               for(j=1; j<(ny[0]-1); j++)
                    for(k=1; k<(nz[0]-1); k++)
                         S[0][i][j][k] = p[0][i][j][k]*epR;
          
          while(t<(niveis-1))
          {
               for(i=1; i<(nx[t+1]-1); i++)
                    for(j=1; j<(ny[t+1]-1); j++)
                         for(k=1; k<(nz[t+1]-1); k++)
                              S[t+1][i][j][k] = 0;
                         
               fronteira(S[t], h[t], nx[t], ny[t], nz[t], D, L, s, fun, fron, 0);                             
               interpola3d(S, nx, ny, nz, t);  
                         
               t++;             
               
               epR = escalar(p[t], R[t], nx[t], ny[t], nz[t]);
                
               for(i=1; i<(nx[t]-1); i++)
                    for(j=1; j<(ny[t]-1); j++) 
                         for(k=1; k<(nz[t]-1); k++)  
                              S[t][i][j][k] += p[t][i][j][k]*epR;                     
          }    
     
          t = niveis-1;    
          corrige(P[t], S[t], nx[t], ny[t], nz[t]);
          fronteira(P[t], h[t], nx[t], ny[t], nz[t], D, L, s, fun, fron, 1); 
          residuo(R[t], P[t], B[t], Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t]);              
          rsd = maximo(R[t], nx[t], ny[t], nz[t]);
          cout << "\nNC = " << NC << " - RESIDUO = " << rsd << endl; 
          
          sor(P[t], B[t], Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t], D, L, s, h[t], w, it, 1, fun, fron);
          
          fprintf(g1, "%i       %e\n", NC, rsd);     
          erro(Pex, P[t], E, nx[t], ny[t], nz[t]);
          emax = maximo(E, nx[t], ny[t], nz[t]);	      
          fprintf(g2, "%i       %e\n", NC, emax);   
             
          if(pre == 's')     
               igual(p[t], R[t], nx[t], ny[t], nz[t]);
          else if(pre == 'j')
               jacobi(p[t], R[t], Ap[t], nx[t], ny[t], nz[t]);            
             
             
          //cin >> pausa;  
             
     }
     while(rsd > eps);
     
     fclose(g1);
     fclose(g2);
     
}


void cg(double ***Pex, double ***P, double ***B, double ***E, double ***z1, double ***z2, double ***R1, double ***R2, double ***p, double ***W, double ****R, double ****rhs, double ****e, double ****Ap, double ****Ae, double ****Aw, double ****An, double ****As, double ****Au, double ****Ad, double ****Mp, double ****Me, double ****Mw, double ****Mn, double ****Ms, double ****Mu, double ****Md, double ****Zero, int *nx, int *ny, int *nz, double D, double L, double S, double *h, double *w, int *itr_d, int *itr_u, double eps, int fun, char fron, char pre, int niveis)
{
     cout << "\n\n                              CG-PRECONDICIONADO\n\n";
          
     int t, i, j, k;
     int NC = 0;
     
     int pausa;
     
     double rsd, alfa, beta, emax;
     
     
     FILE *g1, *g2;
	 g1 = fopen("RSDxNC.dat", "w");	
     g2 = fopen("ERROxNC.dat", "w");
     fprintf(g1, "'NC'       'RSD'\n\n"); 
     fprintf(g2, "'NC'       'EMAX'\n\n");
          
          
     t = niveis-1;
          
     fronteira(P, h[t], nx[t], ny[t], nz[t], D, L, S, fun, fron, 1);      
     residuo(R1, P, B, Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t]);
     rsd = maximo(R1, nx[t], ny[t], nz[t]);     
     cout << "\nNC: " << NC << "     -     RSD: " << rsd << endl;
     
     
     fprintf(g1, "%i       %e\n", NC, rsd);     
     erro(Pex, P, E, nx[t], ny[t], nz[t]);
     emax = maximo(E, nx[t], ny[t], nz[t]);	      
     fprintf(g2, "%i       %e\n", NC, emax);     
     
     
     if(pre == 's')     
          igual(z1, R1, nx[t], ny[t], nz[t]);
     else if(pre == 'j')
          jacobi(z1, R1, Ap[t], nx[t], ny[t], nz[t]);
     else if(pre == 'm')
          MG(Pex, E, niveis, e, R, rhs, z1, R1, Mp, Me, Mw, Mn, Ms, Mu, Md, itr_u, itr_d, nx, ny, nz, D, L, S, h, w, eps, fun, fron);

      
     igual(p, z1, nx[t], ny[t], nz[t]);      
         
     fronteira(p, h[t], nx[t], ny[t], nz[t], D, L, S, fun, fron, 0);   
                   
     for(i=1; i<(nx[t]-1); i++)
          for(j=1; j<(ny[t]-1); j++)
               for(k=1; k<(nz[t]-1); k++)
                    W[i][j][k] = Ap[t][i][j][k]*p[i][j][k] + Ae[t][i][j][k]*p[i+1][j][k] + Aw[t][i][j][k]*p[i-1][j][k] 
                                 + An[t][i][j][k]*p[i][j+1][k] + As[t][i][j][k]*p[i][j-1][k] 
                                 + Au[t][i][j][k]*p[i][j][k+1] + Ad[t][i][j][k]*p[i][j][k-1];
                     
     alfa = escalar(R1, z1, nx[t], ny[t], nz[t])/escalar(p, W, nx[t], ny[t], nz[t]);
    
     do
     {
          NC++;
                    
          for(i=1; i<(nx[t]-1); i++)
               for(j=1; j<(ny[t]-1); j++)
                    for(k=1; k<(nz[t]-1); k++)
                         P[i][j][k] += alfa*p[i][j][k];
                                   
          fronteira(P, h[t], nx[t], ny[t], nz[t], D, L, S, fun, fron, 1);                  
          residuo(R2, P, B, Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t]);
          rsd = maximo(R2, nx[t], ny[t], nz[t]);          
          cout << "\nNC: " << NC << "     -     RSD: " << rsd << endl;
          
          
          fprintf(g1, "%i       %e\n", NC, rsd);     
          erro(Pex, P, E, nx[t], ny[t], nz[t]);
          emax = maximo(E, nx[t], ny[t], nz[t]);	      
          fprintf(g2, "%i       %e\n", NC, emax);
          
          
          if(pre == 's')     
               igual(z2, R2, nx[t], ny[t], nz[t]);
          else if(pre == 'j')
               jacobi(z2, R2, Ap[t], nx[t], ny[t], nz[t]);
          else if(pre == 'm')
               MG(Pex, E, niveis, e, R, rhs, z2, R2, Mp, Me, Mw, Mn, Ms, Mu, Md, itr_u, itr_d, nx, ny, nz, D, L, S, h, w, eps, fun, fron);

          
          beta = escalar(R2, z2, nx[t], ny[t], nz[t])/escalar(R1, z1, nx[t], ny[t], nz[t]);
          
          for(i=1; i<(nx[t]-1); i++)
               for(j=1; j<(ny[t]-1); j++)
                    for(k=1; k<(nz[t]-1); k++)
                         p[i][j][k] = z2[i][j][k] + beta*p[i][j][k];
          
          fronteira(p, h[t], nx[t], ny[t], nz[t], D, L, S, fun, fron, 0);          
                         
          for(i=1; i<(nx[t]-1); i++)
               for(j=1; j<(ny[t]-1); j++)
                    for(k=1; k<(nz[t]-1); k++)
                         W[i][j][k] = Ap[t][i][j][k]*p[i][j][k] + Ae[t][i][j][k]*p[i+1][j][k] + Aw[t][i][j][k]*p[i-1][j][k] 
                                      + An[t][i][j][k]*p[i][j+1][k] + As[t][i][j][k]*p[i][j-1][k] 
                                      + Au[t][i][j][k]*p[i][j][k+1] + Ad[t][i][j][k]*p[i][j][k-1];
                                      
          alfa = escalar(R2, z2, nx[t], ny[t], nz[t])/escalar(p, W, nx[t], ny[t], nz[t]);
          
          igual(R1, R2, nx[t], ny[t], nz[t]);
          igual(z1, z2, nx[t], ny[t], nz[t]); 
                   
          //cin >> pausa;
     }
     while(rsd > eps);
     
     fclose(g1);
     fclose(g2);
     
} 
    
         
void gradiente(double ***Pex, double ***E, double ***P, double ***B, double ***R, double ***W, double ***Ap, double ***Aw, double ***Ae, double ***An, double ***As, double ***Au, double ***Ad, int nx, int ny, int nz, double D, double L, double S, double h, double eps, int fun, char fron)
{
     cout << "\n\n                        GRADIENTE\n\n";
     
     FILE *g1, *g2;
	 g1 = fopen("RSDxNC.dat", "w");	
     g2 = fopen("ERROxNC.dat", "w");
     fprintf(g1, "'NC'       'RSD'\n\n"); 
     fprintf(g2, "'NC'       'EMAX'\n\n"); 
          
     int i, j, k;      
     int pausa;
     int NC = 0; 
                            
     double alfa;
     double rsd, emax;     
                 
     do
     {
          NC++;
                 
          fronteira(P, h, nx, ny, nz, D, L, S, fun, fron, 1);          
          residuo(R, P, B, Ap, Ae, Aw, An, As, Au, Ad, nx, ny, nz);
          rsd = maximo(R, nx, ny, nz);
          cout << "NC: " << NC << " - RSD: " << rsd << "\n";
          fprintf(g1, "%i       %e\n", NC, rsd);
          
          erro(Pex, P, E, nx, ny, nz);
	      emax = maximo(E, nx, ny, nz);	      
	      fprintf(g2, "%i       %e\n", NC, emax);
          

          fronteira(R, h, nx, ny, nz, D, L, S, fun, fron, 0); 
                    
          for(i=1; i<(nx-1); i++)
               for(j=1; j<(ny-1); j++)
                    for(k=1; k<(nz-1); k++)
                         W[i][j][k] = Ap[i][j][k]*R[i][j][k] + Ae[i][j][k]*R[i+1][j][k] + Aw[i][j][k]*R[i-1][j][k] + An[i][j][k]*R[i][j+1][k] + As[i][j][k]*R[i][j-1][k] + Au[i][j][k]*R[i][j][k+1] + Ad[i][j][k]*R[i][j][k-1];
                          
          alfa = escalar(R, R, nx, ny, nz)/escalar(R, W, nx, ny, nz);
                    
          for(i=1; i<(nx-1); i++)
               for(j=1; j<(ny-1); j++)
                    for(k=1; k<(nz-1); k++)
                         P[i][j][k] = P[i][j][k] + alfa*R[i][j][k];
          
          //cin >> pausa;
     }
     while(rsd > eps); 
     fclose(g1);
     fclose(g2);       
}   
*/           
           

                
                
       


     
     
