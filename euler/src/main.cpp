/*

"Busquem, pois, em primeiro lugar o Reino de Deus e a sua justiça, e todas
essas coisas lhes serão acrescentadas." Mt 6.33

*/

//PROJETO 3D NAVIER-STOKES - FUNÇÕES MANUFATURADAS


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
#include "euler.h"
#include "fronteiras.h"

using namespace std;

int main(int argc, char *argv[])
{
    int t, i, j, k;
    
    int pausa;    
    
    int K = 0;
    char fron = 'p';
    
    
    
    //MALHAS (CG): (2 - 4 - 8 - 16 - 32 - 64 - 128 - 256 - 512 - 1024) + 2
    double D = 1;	
    int cgy = 4;
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
    //MALHAS (CG): (2 - 4 - 8 - 16 - 32 - 64 - 128 - 256 - 512 - 1024) + 2
    
    
    
    //DADOS DO MULTIGRID -------------------------------------------------------   	
	const int niveis = 4;

	double *omega = new double [niveis];
	int *itr_u = new int [niveis];
    int *itr_d = new int [niveis];
     
    omega[0] = 1.4;            itr_d[0] = 6;                      itr_u[0] = 6;
    omega[1] = 1.4;            itr_d[1] = 2;                      itr_u[1] = 2;
    omega[2] = 1.4;            itr_d[2] = 2;                      itr_u[2] = 2;
    omega[3] = 1.4;            itr_d[3] = 2;                      itr_u[3] = 2;
    //omega[4] = 1.4;            itr_d[4] = 2;                      itr_u[4] = 2;
    //omega[5] = 1.4;            itr_d[5] = 1;                      itr_u[5] = 1;
    
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
    //DADOS DO MULTIGRID ------------------------------------------------------- 
    
    
    
    t = niveis-1;
    
    
    
    //MOSTRAR PARAMETROS -------------------------------------------------------	 
    double eps = h[t]*h[t];
	 
    cout << "\nNX = " << nx[t]-2 << endl;
    cout << "NY = " << ny[t]-2 << endl;
    cout << "NZ = " << nz[t]-2 << endl;
    cout << "H = " << h[t] << endl;
	 
    cout << "\nEPS = " << eps << "\n\n";
    //MOSTRAR PARAMETROS -------------------------------------------------------
    
    
    
    //PROPRIEDADES DO ESCOAMENTO -----------------------------------------------         
    double Re = 5000;  
    double Ut = 1;
    //double mi = Ut*L/Re;    
    double mi = 1.0;  
    //PROPRIEDADES DO ESCOAMENTO -----------------------------------------------
    
    
    
    //DADOS DA SIMULAÇÃO E PÓS-PROCESSAMENTO -----------------------------------        
    int nt = 2;
    double dt = 1E-03;    
    int NT = ceil(nt/dt);
    int arquivos = 10;       
    //DADOS DA SIMULAÇÃO E PÓS-PROCESSAMENTO -----------------------------------  
  
  
           
    //MATRIZES 3D --------------------------------------------------------------
    double ***Pex = new double **[nx[t]+1];   
    double ***E = new double **[nx[t]+1];     
    double ***P = new double **[nx[t]+1];   
     
    double ***Fx = new double **[nx[t]+1];
    double ***Fy = new double **[nx[t]+1];
    double ***Fz = new double **[nx[t]+1];    
    
    double ***Uex = new double **[nx[t]+1]; 
    double ***Ue = new double **[nx[t]+1]; 
    double ***U = new double **[nx[t]+1]; 
    
    double ***Vex = new double **[nx[t]+1]; 
    double ***Ve = new double **[nx[t]+1]; 
    double ***V = new double **[nx[t]+1]; 
    
    double ***Wex = new double **[nx[t]+1]; 
    double ***We = new double **[nx[t]+1]; 
    double ***W = new double **[nx[t]+1]; 
    
    for(i=0; i<(nx[t]+1); i++)
    {
         Pex[i] = new double *[ny[t]+1];  
         E[i] = new double *[ny[t]+1];            
         P[i] = new double *[ny[t]+1];  
                    
         Fx[i] = new double *[ny[t]+1];
         Fy[i] = new double *[ny[t]+1];
         Fz[i] = new double *[ny[t]+1];  
         
         Uex[i] = new double *[ny[t]+1]; 
         Ue[i] = new double *[ny[t]+1]; 
         U[i] = new double *[ny[t]+1]; 
         
         Vex[i] = new double *[ny[t]+1]; 
         Ve[i] = new double *[ny[t]+1]; 
         V[i] = new double *[ny[t]+1]; 
         
         Wex[i] = new double *[ny[t]+1]; 
         We[i] = new double *[ny[t]+1]; 
         W[i] = new double *[ny[t]+1];  
                   
         for(j=0; j<(ny[t]+1); j++)
         {
              Pex[i][j] = new double [nz[t]+1]; 
              E[i][j] = new double [nz[t]+1];             
              P[i][j] = new double [nz[t]+1];  
                          
              Fx[i][j] = new double [nz[t]+1];
              Fy[i][j] = new double [nz[t]+1];
              Fz[i][j] = new double [nz[t]+1];  
              
              Uex[i][j] = new double [nz[t]+1];  
              Ue[i][j] = new double [nz[t]+1];  
              U[i][j] = new double [nz[t]+1];  
              
              Vex[i][j] = new double [nz[t]+1];  
              Ve[i][j] = new double [nz[t]+1];  
              V[i][j] = new double [nz[t]+1];  
              
              Wex[i][j] = new double [nz[t]+1];  
              We[i][j] = new double [nz[t]+1];  
              W[i][j] = new double [nz[t]+1];  
                          
              for(k=0; k<(nz[t]+1); k++)
              {
                   Pex[i][j][k] = 0;  
                   E[i][j][k] = 0;                
                   P[i][j][k] = 0; 
                                 
                   Fx[i][j][k] = 0;
                   Fy[i][j][k] = 0;
                   Fz[i][j][k] = 0;
                   
                   Uex[i][j][k] = 0;
                   Ue[i][j][k] = 0;
                   U[i][j][k] = 0;
                   
                   Vex[i][j][k] = 0;
                   Ve[i][j][k] = 0;
                   V[i][j][k] = 0;
                   
                   Wex[i][j][k] = 0;
                   We[i][j][k] = 0;
                   W[i][j][k] = 0;
              }
         }
    }   
    //MATRIZES 3D --------------------------------------------------------------  
    
        
       
    //MATRIZES 4D --------------------------------------------------------------
	double ****Ap = new double ***[niveis];
    double ****Ae = new double ***[niveis];
    double ****Aw = new double ***[niveis];
    double ****An = new double ***[niveis];
    double ****As = new double ***[niveis];
	double ****Au = new double ***[niveis];
	double ****Ad = new double ***[niveis];
	double ****Ro = new double ***[niveis];
		
	double ****Pe = new double ***[niveis];	
    double ****B = new double ***[niveis];
    
    double ****R = new double ***[niveis];	 
    double ****e = new double ***[niveis];     
    double ****rhs = new double ***[niveis];
    
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
			
		Pe[t] = new double **[nx[t]];	
        B[t] = new double **[nx[t]];
          
        R[t] = new double **[nx[t]];          
        e[t] = new double **[nx[t]];
        rhs[t] = new double **[nx[t]];
          
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
	
			Pe[t][i] = new double *[ny[t]];
			B[t][i] = new double *[ny[t]];
			   
			R[t][i] = new double *[ny[t]];			   
			e[t][i] = new double *[ny[t]];
			rhs[t][i] = new double *[ny[t]];
			   
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
											    
				Pe[t][i][j] = new double [nz[t]];			
				B[t][i][j] = new double [nz[t]];
				    
				R[t][i][j] = new double [nz[t]];				    
				e[t][i][j] = new double [nz[t]];
				rhs[t][i][j] = new double [nz[t]];
				    
				for(k=0; k<nz[t]; k++)
				{
                    Ap[t][i][j][k] = 0;
                    Ae[t][i][j][k] = 0;
       	            Aw[t][i][j][k] = 0;
                    An[t][i][j][k] = 0;
                    As[t][i][j][k] = 0;
					Au[t][i][j][k] = 0;
					Ad[t][i][j][k] = 0;
					Ro[t][i][j][k] = 0;
									     
					Pe[t][i][j][k] = 0;			
					B[t][i][j][k] = 0;
					     
					R[t][i][j][k] = 0;					     
					e[t][i][j][k] = 0;
					rhs[t][i][j][k] = 0;
                }
            }
        }
    }
    //MATRIZES 4D --------------------------------------------------------------
    
    
    
    ME1(Ro, nx, ny, nz, h, K, 'n', niveis);
    coeficientes(Ap, Aw, Ae, An, As, Au, Ad, Ro, h, nx, ny, nz, niveis);
    
    
    
    t = niveis-1;
    
    
    
    //TEMPO INICIAL 0
    P_EX(P, h[t], nx[t], ny[t], nz[t], 'n', 0.0);
    U_EX(U, h[t], nx[t], ny[t], nz[t], 'n', 0.0);
    V_EX(V, h[t], nx[t], ny[t], nz[t], 'n', 0.0);  
    W_EX(W, h[t], nx[t], ny[t], nz[t], 'n', 0.0);
    
    
        
    double rsd, emax;
    int tt = 0;
    
    FILE *comp;
    comp = fopen("NCxTempo.dat","w");
    fprintf(comp, "Tempo          NC\n\n");
    
    cin >> pausa;
     
     
    for(int tempo=1; tempo<=10000; tempo++)
    {       
         
         FX(Fx, Ro[t], mi, h[t], nx[t], ny[t], nz[t], 'n', tempo*dt);
         FY(Fy, Ro[t], mi, h[t], nx[t], ny[t], nz[t], 'n', tempo*dt);
         FZ(Fz, Ro[t], mi, h[t], nx[t], ny[t], nz[t], 'n', tempo*dt);    
         
         
         P_EX(Pex, h[t], nx[t], ny[t], nz[t], 'n', tempo*dt);
         U_EX(Uex, h[t], nx[t], ny[t], nz[t], 'n', tempo*dt);
         V_EX(Vex, h[t], nx[t], ny[t], nz[t], 'n', tempo*dt);  
         W_EX(Wex, h[t], nx[t], ny[t], nz[t], 'n', tempo*dt);                
        
         //Dirichlet(nx[t], ny[t], nz[t], U, V, W, Uex, Vex, Wex);
         Fronteira_P(P, h[t], nx[t], ny[t], nz[t], 'p', 1);
         Fronteira_U(U, h[t], nx[t], ny[t], nz[t], 'p');
         Fronteira_V(V, h[t], nx[t], ny[t], nz[t], 'p');
         Fronteira_W(W, h[t], nx[t], ny[t], nz[t], 'p');
         
         
         //Cavidade_vel(nx[t], ny[t], nz[t], U, V, W, Ut);



         ue_laminar(nx[t], ny[t], nz[t], dt, Ro[t], mi, Ue, U, V, W, P, Fx, h[t]);  
         ve_laminar(nx[t], ny[t], nz[t], dt, Ro[t], mi, Ve, U, V, W, P, Fy, h[t]);
         we_laminar(nx[t], ny[t], nz[t], dt, Ro[t], mi, We, U, V, W, P, Fz, h[t]);
   


         Update_vel(nx[t], ny[t], nz[t], Ue, Ve, We, Uex, Vex, Wex);
         //Update_vel_cavidade(nx[t], ny[t], nz[t], Ue, Ve, We);



         f(nx[t], ny[t], nz[t], B[t], Ue, Ve, We, dt, h[t]); 
         normaliza(B[t], nx[t], ny[t], nz[t], h[t]);     
         //P_EX(P, h[t], nx[t], ny[t], nz[t], 'n', tempo*dt);
         //P_EX(Pe[t], h[t], nx[t], ny[t], nz[t], 'n', tempo*dt);
         
         
         tt = MG(niveis, e, R, rhs, Pe[t], B[t], Ap, Ae, Aw, An, As, Au, Ad, itr_u, itr_d, nx, ny, nz, D, L, S, h, omega, eps, fron);


         
         /*         
         do
         {
              tt++;
              sor(Pe[t], B[t], Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t], D, L, S, h[t], omega[t], 1, 1, fron);
              residuo(R[t], Pe[t], B[t], Ap[t], Ae[t], Aw[t], An[t], As[t], Au[t], Ad[t], nx[t], ny[t], nz[t]);
              rsd = maximo(R[t], nx[t], ny[t], nz[t]);
                           
              cout << "\nNC: " << tt << "      -      rsd: " << rsd;                           
         }
         while(rsd > eps);
         */
         
         
         
         fprintf(comp, "%i            %i\n", tempo, tt);
         tt = 0;  
         
         cout << "\n\n";                
         
         
                        
         p(nx[t], ny[t], nz[t], P, Pe[t]);         
         u(nx[t], ny[t], nz[t], dt, Ro[t], U, Ue, Pe[t], h[t]);         
         v(nx[t], ny[t], nz[t], dt, Ro[t], V, Ve, Pe[t], h[t]);         
         w(nx[t], ny[t], nz[t], dt, Ro[t], W, We, Pe[t], h[t]);
         
         
         //zero(Pe, nx, ny, nz, niveis);
         
         continuidade(U, V, W, nx[t], ny[t], nz[t], h[t], tempo);

                
                        
         if( fmod(tempo, ceil(1/(dt*arquivos))) == 0 )           
              campo(nx[t], ny[t], nz[t], h[t], U, V, W, P, tempo);               
                     
                     
             
         //ERRO    
         erro(Pex, P, E, nx[t], ny[t], nz[t]);
	     emax = maximo(E, nx[t], ny[t], nz[t]);
	     cout << "\n\nemax_P = " << emax << endl;
	
	     erro(Uex, U, E, nx[t], ny[t], nz[t]);
	     emax = maximo(E, nx[t], ny[t], nz[t]);
	     cout << "\nemax_U = " << emax << endl;
	
	     erro(Vex, V, E, nx[t], ny[t], nz[t]);
	     emax = maximo(E, nx[t], ny[t], nz[t]);
	     cout << "\nemax_V = " << emax << endl;
	
	     erro(Wex, W, E, nx[t], ny[t], nz[t]);
	     emax = maximo(E, nx[t], ny[t], nz[t]);
	     cout << "\nemax_W = " << emax << endl;            
                     
                   
                     
                      
        // cin >> pausa;
                            
    }
    fclose(comp);
    
  
    
    
    
    system("PAUSE");
    return EXIT_SUCCESS;
}
