/*



EQUAÇÃO DA ENERGIA BI-DIMENSIONAL



*/



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
#include "impressao.h"


using namespace std;

int main(int argc, char *argv[])
{    
     cout << "ENERGIA...\n\n\n";     
     
     int pausa;              
     int t, i, j;             
        
	 //MALHA GROSSEIRA: (2 - 4 - 8 - 16 - 32 - 64 - 128 - 256 - 512 - 1024) + 2
     	 
	 double D = 1;	
	 int cgy = 4;
	 double H = D/(cgy-2);
	 int mult = 1;
	 int cgx = cgy-2;
	 for(int m=1; m<mult; m++)	 
          cgx = 2*cgx;          
     cgx = cgx+2;
	 double L = pow(2.0,(mult-1))*D;
	 
	 cout << "D = " << D << "     -     (CGY-2) = " << cgy-2 << endl;
	 cout << "\nMULT = " << mult << endl;
	 cout << "\nL = " << L << "     -     (CGX-2) = " << cgx-2 << endl;
	 		
	 const int niveis = 5;
	 
	 int temp_nx = cgx-2;
	 int temp_ny = cgy-2;

	 int *nx = new int [niveis];
	 int *ny = new int [niveis];
	 double *h = new double [niveis];

	 for(t=0; t<niveis; t++)
	 {
          h[t] = D/temp_ny;
		  nx[t] = temp_nx+2;
		  ny[t] = temp_ny+2;
		  temp_nx = 2*temp_nx;
		  temp_ny = 2*temp_ny;		
	 } 
	 
	 
	 double *w = new double [niveis];    int *itr_u = new int [niveis];    int *itr_d = new int [niveis];
     w[0] = 1;                           itr_u[0] = 10;                    itr_d[0] = 10;
     w[1] = 1;                           itr_u[1] = 2;                     itr_d[1] = 2;
     w[2] = 1;                           itr_u[2] = 2;                     itr_d[2] = 2;
     w[3] = 1;                           itr_u[3] = 2;                     itr_d[3] = 2;
     w[4] = 1;                           itr_u[4] = 2;                     itr_d[4] = 2;
     //w[5] = 1;                           itr_u[5] = 2;                     itr_d[5] = 2; 
     //w[6] = 1;                           itr_u[6] = 2;                     itr_d[6] = 2;
     
         
	 
	 t = niveis-1;
	 
	 	 
	 double eps = h[t]*h[t];
     double dt = 1E-03; 
           
     cout << "\nNIVEIS = " << niveis << endl;	 
	 cout << "\nNX = " << nx[t]-2 << endl;  
	 cout << "\nNY = " << ny[t]-2 << endl;
	 cout << "\nh = " << h[t] << endl;
     cout << "\nEPS = " << eps << endl; 
     cout << "\ndt = " << dt << "\n\n";
     
     
     cout << "\n\nPAUSA!" << endl;
     cin >> pausa;
          
     
     //DADOS DO PROBLEMA     
     double Tw = 1400;
     double Te = 1000;
     double Ti = 800;
     double ro = 7200;
     double Latente = -265200;
     
     int arquivos = 2; 
        
     
     //MATRIZES   
     double ***Ap = new double **[niveis];
     double ***Ae = new double **[niveis];
     double ***Aw = new double **[niveis];
     double ***An = new double **[niveis];
     double ***As = new double **[niveis];
     double ***Cp = new double **[niveis];
     double ***K = new double **[niveis];      
     
     double ***e = new double **[niveis];
     double ***R = new double **[niveis];
     double ***rhs = new double **[niveis];
          
     double ***R1 = new double **[niveis];
     double ***R2 = new double **[niveis];
     double ***z1 = new double **[niveis];
     double ***z2 = new double **[niveis];
     double ***p = new double **[niveis];
     double ***W = new double **[niveis];    
     
     for(t=0; t<niveis; t++)
     {
	      Ap[t] = new double *[nx[t]];
          Ae[t] = new double *[nx[t]];
          Aw[t] = new double *[nx[t]];
          An[t] = new double *[nx[t]];
          As[t] = new double *[nx[t]];         
          Cp[t] = new double *[nx[t]];          
          K[t] = new double *[nx[t]];          
          
          e[t] = new double *[nx[t]];
          R[t] = new double *[nx[t]];
          rhs[t] = new double *[nx[t]];      
                    
          R1[t] = new double *[nx[t]];
          R2[t] = new double *[nx[t]];
          z1[t] = new double *[nx[t]];
          z2[t] = new double *[nx[t]];
          p[t] = new double *[nx[t]];
          W[t] = new double *[nx[t]];
                
          for(i=0; i<nx[t]; i++)
          {
               Ap[t][i] = new double [ny[t]];
               Ae[t][i] = new double [ny[t]];
               Aw[t][i] = new double [ny[t]];
               An[t][i] = new double [ny[t]];
               As[t][i] = new double [ny[t]];               
               Cp[t][i] = new double [ny[t]];               
               K[t][i] = new double [ny[t]];               
               
               e[t][i] = new double [ny[t]];
               R[t][i] = new double [ny[t]];
               rhs[t][i] = new double [ny[t]];          
               
               R1[t][i] = new double [ny[t]];
               R2[t][i] = new double [ny[t]];
               z1[t][i] = new double [ny[t]];
               z2[t][i] = new double [ny[t]];
               p[t][i] = new double [ny[t]];
               W[t][i] = new double [ny[t]];       
               
               for(j=0; j<ny[t]; j++)
               {
		            Ap[t][i][j] = 0;
                    Ae[t][i][j] = 0;
                    Aw[t][i][j] = 0;
				    An[t][i][j] = 0;
                    As[t][i][j] = 0;                    
                    Cp[t][i][j] = 0;                    
                    K[t][i][j] = 0;                   
                    
                    e[t][i][j] = 0;
                    R[t][i][j] = 0;
                    rhs[t][i][j] = 0;                    
                    
                    R1[t][i][j] = 0;
                    R2[t][i][j] = 0;
                    z1[t][i][j] = 0;
                    z2[t][i][j] = 0;
                    p[t][i][j] = 0;
                    W[t][i][j] = 0;                   
			   }
          }
     }
     
     t = niveis-1;
     
     double **T = new double *[nx[t]];
     double **Told = new double *[nx[t]];     
     double **Cp0 = new double *[nx[t]];
     double **f = new double *[nx[t]];
     double **f0 = new double *[nx[t]];
     double **B = new double *[nx[t]]; 
     char **flag = new char *[nx[t]];
     
     for(i=0; i<nx[t]; i++)
     {
          T[i] = new double [ny[t]]; 
          Told[i] = new double [ny[t]];
          Cp0[i] = new double [ny[t]];
          f[i] = new double [ny[t]];
          f0[i] = new double [ny[t]];
          B[i] = new double [ny[t]];
          flag[i] = new char [ny[t]];
          
          for(j=0; j<ny[t]; j++)
          {
               T[i][j] = 0;   
               Told[i][j] = 0;
               Cp0[i][j] = 0;
               f[i][j] = 0;
               f0[i][j] = 0;
               B[i][j] = 0;
          }
     }
     
     Set_Flag(flag, nx[t], ny[t]);
         
     Set_Initial_Temperature(T, Told, nx[t], ny[t], Ti);    
                  
     int totalTime = 20;
     int NC = 0;
     double rsd;
          
     cout << "\n\nPASSOS TOTAIS = " << totalTime/dt << "\n\n";
                   
     for(int time=0; time<=ceil(totalTime/dt); time++)
     {    
          cout << "time = " << time << endl;
                    
          set_K(K, T, nx, ny, niveis);
          
          set_Cp0(Cp0, T, nx[t], ny[t]);
          
          set_Cp(Cp, Cp0, T, Told, nx, ny, niveis);
          
          set_f0(f0, T, nx[t], ny[t]);
          
          set_f(f, f0, T, Told, nx[t], ny[t]);
          
          igual(Told, T, nx[t], ny[t]);
           
          coefficients(flag, Ap, Ae, Aw, An, As, B, Cp, K, T, Cp0, f, f0, nx, ny, h, ro, Latente, dt, niveis, Tw, Te); 
 
          t = niveis-1;
          
          if( (time == 0) || (fmod(time, ceil(1/(dt*arquivos))) == 0) )            
               campo(nx[t], ny[t], h[t], Cp[t], K[t], T, f, time);
          
          do
          {               
               NC++;
               sor(T, K[t], B, Ap[t], Ae[t], Aw[t], An[t], As[t], nx[t], ny[t], L, D, h[t], w[t], 1);
               residuo(R[t], T, B, Ap[t], Ae[t], Aw[t], An[t], As[t], nx[t], ny[t]);
               rsd = maximo(R[t], nx[t], ny[t]);
               cout << "NC: " << NC << "      -      rsd: " << rsd << endl;               
          }
          while(rsd > eps);
          
          NC = 0;
          
          //cin >> pausa;  
              
          
     }   
     
     
          
    
     system("PAUSE");
     return EXIT_SUCCESS;
}
