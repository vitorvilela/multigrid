//DOMINIO QUADRADO

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
#include "funcoes.h"

using namespace std;

int main(int argc, char *argv[])
{
    
     cout << "\n\n                       INICIO DA SIMULACAO\n\n";
    
     int pausa;       
    
     int t, i, j, k;   
      
    
     const int fun = 1;
     const char fron = 'n';
     const char pre = 'j';
     
     const int K = 1000;  
         
    
	 //cgy = (2 - 4 - 8 - 16 - 32 - 64 - 128 - 256) + 2
	 double D = 1;	
	 int cgy = 6;
	 double H = D/(cgy-2);
	 
	 int mult_x = 1;
	 int mult_z = 1;
	 
     int cgx = cgy-2;
     int cgz = cgy-2;
     
	 for(int m=1; m<mult_x; m++)
          cgx = 2*cgx;
     cgx = cgx+2;
     for(int m=1; m<mult_z; m++)
          cgz = 2*cgz;
     cgz = cgz+2;          
          
	 double L = pow((cgy-2),(mult_x-1))*D;	
	 double S = pow((cgy-2),(mult_z-1))*D;
	
	 
	 //DADOS DO MULTIGRID     	
	 const int niveis = 4;

	 double *w = new double [niveis];
	 int *itr_u = new int [niveis];
     int *itr_d = new int [niveis];
     
     w[0] = 1;            itr_d[0] = 6;                      itr_u[0] = 6;
     w[1] = 1;            itr_d[1] = 2;                       itr_u[1] = 2;
     w[2] = 1;            itr_d[2] = 2;                       itr_u[2] = 2;
     w[3] = 1;            itr_d[3] = 2;                       itr_u[3] = 2;
     //w[4] = 1;            itr_d[4] = 2;                       itr_u[4] = 2;
     //w[5] = 1;            itr_d[5] = 1;                       itr_u[5] = 1;
     //w[6] = 1;            itr_d[6] = 1;                       itr_u[6] = 1;
     //w[7] = 1;            itr_d[7] = 1;                       itr_u[7] = 1;
     //w[8] = 1;            itr_d[8] = 1;                       itr_u[8] = 1;
     //w[9] = 1;            itr_d[9] = 1;                       itr_u[9] = 1;
    
                   
    
	 //SUBNIVEIS
	 int temp_nx = cgx-2;
	 int temp_ny = cgy-2;
	 int temp_nz = cgz-2;

	 int *nx = new int [niveis];
	 int *ny = new int [niveis];
	 int *nz = new int [niveis];
	 double *h = new double [niveis];

	 for(t=0; t<niveis; t++)
	 {
          h[t] = D/temp_ny;
		  nx[t] = temp_nx+2;
		  ny[t] = temp_ny+2;
		  nz[t] = temp_nz+2;
		  temp_nx = 2*temp_nx;
		  temp_ny = 2*temp_ny;
          temp_nz = 2*temp_nz;		
	 }
	 
	 
	 	 
	 //MATRIZES 4D
	 double ****Ap = new double ***[niveis];
     double ****Ae = new double ***[niveis];
     double ****Aw = new double ***[niveis];
     double ****An = new double ***[niveis];
     double ****As = new double ***[niveis];
	 double ****Au = new double ***[niveis];
	 double ****Ad = new double ***[niveis];
	 double ****Ro = new double ***[niveis];	 
	 
	 double ****Pex = new double ***[niveis];
	 double ****P = new double ***[niveis];
     double ****B = new double ***[niveis];
	 double ****E = new double ***[niveis];
     	
	 double ****R = new double ***[niveis];	 
	 double ****e = new double ***[niveis];     
     double ****rhs = new double ***[niveis]; 
     
     double ****z1 = new double ***[niveis];
     double ****z2 = new double ***[niveis];
     double ****R1 = new double ***[niveis];
     double ****R2 = new double ***[niveis];         
     double ****p = new double ***[niveis];
     double ****W = new double ***[niveis];
     
     double ****Zero = new double ***[niveis];
     	 
     for(t=0; t<niveis; t++)
     {
          Ap[t] = new double **[nx[t]];
          Ae[t] = new double **[nx[t]];
          Aw[t] = new double **[nx[t]];
          An[t] = new double **[nx[t]];
          As[t] = new double **[nx[t]];
		  Au[t] = new double **[nx[t]];
		  Ad[t] = new double **[nx[t]];
		  Ro[t] = new double **[nx[t]];		  
		  
		  Pex[t] = new double **[nx[t]];
          P[t] = new double **[nx[t]];
          B[t] = new double **[nx[t]];
          E[t] = new double **[nx[t]];
          
          R[t] = new double **[nx[t]];          
          e[t] = new double **[nx[t]];
          rhs[t] = new double **[nx[t]];
          
          z1[t] = new double **[nx[t]];
          z2[t] = new double **[nx[t]];
          R1[t] = new double **[nx[t]];
          R2[t] = new double **[nx[t]];
          p[t] = new double **[nx[t]];          
          W[t] = new double **[nx[t]];	 
          
          Zero[t] = new double **[nx[t]]; 
		  
          for(i=0; i<nx[t]; i++)
          {
               Ap[t][i] = new double *[ny[t]];
               Ae[t][i] = new double *[ny[t]];
               Aw[t][i] = new double *[ny[t]];
               An[t][i] = new double *[ny[t]];
               As[t][i] = new double *[ny[t]];
			   Au[t][i] = new double *[ny[t]];
			   Ad[t][i] = new double *[ny[t]];
			   Ro[t][i] = new double *[ny[t]];			   
			   
			   Pex[t][i] = new double *[ny[t]];
			   P[t][i] = new double *[ny[t]];
			   B[t][i] = new double *[ny[t]];
			   E[t][i] = new double *[ny[t]];
			   
			   R[t][i] = new double *[ny[t]];			   
			   e[t][i] = new double *[ny[t]];
			   rhs[t][i] = new double *[ny[t]];
			   
			   z1[t][i] = new double *[ny[t]];
			   z2[t][i] = new double *[ny[t]];
			   R1[t][i] = new double *[ny[t]];
			   R2[t][i] = new double *[ny[t]];
			   p[t][i] = new double *[ny[t]];
			   W[t][i] = new double *[ny[t]];
			   
			   Zero[t][i] = new double *[ny[t]];
			   
               for(j=0; j<ny[t]; j++)
               {
				    Ap[t][i][j] = new double [nz[t]];
				    Ae[t][i][j] = new double [nz[t]];
				    Aw[t][i][j] = new double [nz[t]];
				    An[t][i][j] = new double [nz[t]];
				    As[t][i][j] = new double [nz[t]];
				    Au[t][i][j] = new double [nz[t]];
				    Ad[t][i][j] = new double [nz[t]];
				    Ro[t][i][j] = new double [nz[t]];
				   				    
				    Pex[t][i][j] = new double [nz[t]];
				    P[t][i][j] = new double [nz[t]];
				    B[t][i][j] = new double [nz[t]];
				    E[t][i][j] = new double [nz[t]];
				    
				    R[t][i][j] = new double [nz[t]];				    
				    e[t][i][j] = new double [nz[t]];
				    rhs[t][i][j] = new double [nz[t]];
				    
				    z1[t][i][j] = new double [nz[t]];
				    z2[t][i][j] = new double [nz[t]];
				    R1[t][i][j] = new double [nz[t]];
				    R2[t][i][j] = new double [nz[t]];
				    p[t][i][j] = new double [nz[t]];
				    W[t][i][j] = new double [nz[t]];
				    
				    Zero[t][i][j] = new double [nz[t]];
				    
				    for(k=0; k<nz[t]; k++)
				    {
                         Ap[t][i][j][k] = 0;
                         Ae[t][i][j][k] = 0;
       	                 Aw[t][i][j][k] = 0;
                         An[t][i][j][k] = 0;
                         As[t][i][j][k] = 0;
					     Au[t][i][j][k] = 0;
					     Ad[t][i][j][k] = 0;
					     Ro[t][i][j][k] = 1;
										     
					     Pex[t][i][j][k] = 0;
					     P[t][i][j][k] = 0;
					     B[t][i][j][k] = 0;
					     E[t][i][j][k] = 0;
					     
					     R[t][i][j][k] = 0;					     
					     e[t][i][j][k] = 0;
					     rhs[t][i][j][k] = 0;
					     
                         z1[t][i][j][k] = 0;	
                         z2[t][i][j][k] = 0;
                         R1[t][i][j][k] = 0;
                         R2[t][i][j][k] = 0;
                         p[t][i][j][k] = 0;
                         W[t][i][j][k] = 0;	
                         
                         Zero[t][i][j][k] = 0;				     
				    }
               }
          }
     } 
        
	 ME(Ro, nx, ny, nz, h, K, 'n', niveis);
     coeficientes(Ap, Aw, Ae, An, As, Au, Ad, Ro, h, nx, ny, nz, niveis);
      
     for(t=0; t<niveis; t++)   
     {
	      E1(Pex[t], h[t], nx[t], ny[t], nz[t], 's');
	      f(B[t], Pex[t], Ap[t], Aw[t], Ae[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t]);
     }
    
        
     t = niveis-1;
     cout << "\nMALHA = " << nx[t]-2 << endl;     
     double eps = 1E-08;     
     cout << "\nEPS = " << eps << endl;
     cin >> pausa; 
	 
	 
	 //SOLVER	  
	 time_t start, end;
     double dif;
     time (&start);
	 
	 
	 
     nested(Pex, E, e, R, rhs, z1, z2, R1, R2, p, W, Ap, Ae, Aw, As, An, Au, Ad, Zero, Zero, Zero, Zero, Zero, Zero, Zero, B, P, nx, ny, nz, h, niveis, fun, fron, L, D, S, itr_u, itr_d, eps, w);
	 
	 //MG(Pex[t], E[t], niveis, e, R, rhs, P[t], B[t], Ap, Ae, Aw, An, As, Au, Ad, itr_u, itr_d, nx, ny, nz, D, L, S, h, w, eps, fun, fron);

	 //MGgM(Pex[t], E[t], p, R, rhs, e, Ap, Ae, Aw, As, An, Au, Ad, B, P, nx, ny, nz, h, niveis, fun, fron, L, D, S, eps, pre);

     //cg(Pex[t], P[t], B[t], E[t], z1[t], z2[t], R1[t], R2[t], p[t], W[t], R, rhs, e, Ap, Ae, Aw, An, As, Au, Ad, Zero, Zero, Zero, Zero, Zero, Zero, Zero, Zero, nx, ny, nz, D, L, S, h, w, itr_d, itr_u, eps, fun, fron, pre, niveis);

     //gradiente(Pex[t], E[t], P[t], B[t], R[t], p[t], Ap[t], Aw[t], Ae[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t], D, L, S, h[t], eps, fun, fron);


/*
	 int tempo = 0;	 
	 double rsd = 1;
     double emax = 1;
	 
	 FILE *g1, *g2;
	 g1 = fopen("RSDxNC.dat", "w");	
     g2 = fopen("ERROxNC.dat", "w");
     fprintf(g1, "'NC'       'RSD'\n\n"); 
     fprintf(g2, "'NC'       'EMAX'\n\n"); 
     
     do
     {
          tempo++;
          sor(P[t], B[t], Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t], D, L, S, h[t], w[t], 1, 1, fun, fron);          
          residuo(R[t], P[t], B[t], Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t]);
          rsd = maximo(R[t], nx[t], ny[t], nz[t]);
          
          fprintf(g1, "%i       %e\n", tempo, rsd);          
          cout << "\nNC: " << tempo << "      -      rsd: " << rsd << "\n";          
          
          erro(Pex[t], P[t], E[t], nx[t], ny[t], nz[t]);
	      emax = maximo(E[t], nx[t], ny[t], nz[t]);          	      
	      fprintf(g2, "%i       %e\n", tempo, emax);         
     }
     while(rsd > eps);
      
     fclose(g1);
     fclose(g2);
*/     
    

     time (&end);
     dif = difftime (end, start);
     cout << "\n\nTEMPO = " << dif << "\n";

   
     //ERRO
     double emax;
     erro(Pex[t], P[t], E[t], nx[t], ny[t], nz[t]);
	 emax = maximo(E[t], nx[t], ny[t], nz[t]);
	 cout << "\n\nemax = " << emax << endl;
	 
	  
    
    
     system("PAUSE");
     return EXIT_SUCCESS;
     
}




