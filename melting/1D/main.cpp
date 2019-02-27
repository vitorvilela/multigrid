/*

UNIDIMENSIONAL AND TRANSIENT TEMPERATURE FIELD SIMULATION

BY: VITOR MACIEL VILELA FERREIRA - UFU

"SE É SOMENTE PARA ESTA VIDA QUE TEMOS ESPERANÇA EM CRISTO,
SOMOS, DE TODOS OS HOMENS, OS MAIS DIGNOS DE COMPAIXÃO." 1Co 15.19

*/

#include <cstdlib>
#include <iostream>
using std::cout;
using std::cin;
using std::endl;
#include <math.h>
#include <stdio.h>

using namespace std;

int main(int argc, char *argv[])
{
    cout << "INITIATING THE PROGRAM" << endl;
    
    int pause, i, time;
       
    //DOMAIN'S CONSTANTS
    const double D = 0.5;
    int n = 64;
    const int ghost = 2;
    n += ghost;
    const double h = D/(n-ghost);
    
    cout << endl;
    cout << "DOMAIN LENGHT: " << D << endl;
    cout << "GRID: " << n-ghost << endl;
    cout << "MESH: " << h << endl;
    
    //VECTORS
    double *T = new double [n];
    double *Told = new double [n];
    double *Ap = new double [n];
    double *Aw = new double [n];
    double *Ae = new double [n];
    double *B = new double [n];
    double *R = new double [n];
    double *K = new double [n];
    double *f = new double [n];
    double *f0 = new double [n];
    double *Cp = new double [n];
    double *Cp0 = new double [n];
    for(i=0; i<n; i++)
    {
        T[i] = 0;
        Told[i] = 0;
        Ap[i] = 0;
        Aw[i] = 0;
        Ae[i] = 0;
        B[i] = 0;
        R[i] = 0;
        K[i] = 0;
        f[i] = 0;
        f0[i] = 0;
        Cp[i] = 0;
        Cp0[i] = 0;
    }
    
    //PROBLEM'S CONSTANTS
    const double L = -265200.0;
    const double ro = 7200.0;
    const double Ts = 1400.0;
    const double Tl = 1455.0;
    const double Tb = 1450.0;    
    const double Ti = 1200.0;    
    double Kw, Ke;
    
    cout << endl;
    cout << "LATENT HEAT: " << L << endl;
    cout << "SPECIFIC MASS: " << ro << endl;
    cout << "SOLID TEMPERATURE: " << Ts << endl;
    cout << "LIQUID TEMPERATURE: " << Tl << endl;
    cout << "BOUNDARY TEMPERATURE: " << Tb << endl;
    cout << "INITIAL TEMPERATURE: " << Ti << endl;
    
    //SYSTEM CONSTANTS
    const int totalTime = 1600;
    const double dt = 0.01;
    const double limitTime = ceil(totalTime/dt);
    const double omega = 1.0;    
    double maximumResidue = 0.0;
    const double eps = 1E-03;
    double iteration = 0.0;
    
    //PRINTING CONSTANTS
    int data = 20;
    int j; 
    double vadx, vady;
    const int nx = n - ghost; 
    const int ny = (n - ghost)/4; 
    
    cout << endl;
    cout << "TOTAL TIME OF SIMULATION: " << totalTime << endl;
    cout << "TIME STEP: " << dt << endl;
    cout << "COEFFICIENT OF S.O.R: " << omega << endl;
    cout << "EPS: " << eps << endl;
    
    //INITIAL TEMPERATURE FIELD
    for(i=0; i<n; i++)
    {
        T[i] = Ti;
        Told[i] = Ti;
    }
    
    cin >> pause;
    
    FILE *sai;
    sai = fopen("PERFIL.dat","w");
    fprintf(sai, "TEMPO             TEMPERATURA            F \n\n");
    
    //LOOP OF TIME
    for(time=0; time<=limitTime; time++)
    {
        cout << endl;
        cout << "TIME: " << time << endl;                
                
        //THERMAL CONDUCTIVITY
        for(i=1; i<(n-1); i++)
        {
            if(T[i] <= Ts)
                K[i] = 14.3 + 0.01983*T[i] - 5.451E-06*pow(T[i],2.0);
            else
                K[i] = 31.37804;
        }        
        K[0] = K[1];
        K[n-1] = K[n-2];                    
        
        //SPECIFIC HEAT - CURRENT TIME
        for(i=1; i<(n-1); i++)
        {
            if(T[i] <= Ts)
                Cp0[i] = 460.5 + 0.4257*T[i] - 5.05E-04*pow(T[i],2.0) + 2.6608E-07*pow(T[i],3.0);
            else
                Cp0[i] = 796.584;
        }
                
        //SPECIFIC HEAT - FUTURE TIME
        for(i=1; i<(n-1); i++)
        {
            if(T[i] <= Ts)
                Cp[i] = Cp0[i] + (0.4257 - 1.01E-03*Told[i] + 7.9824E-07*pow(Told[i],2.0))*(T[i]-Told[i]);
            else
                Cp[i] = Cp0[i];
        }
                
        //LIQUID MASS FRACTION - CURRENT TIME
        for(i=1; i<(n-1); i++)
        {
            if(T[i] < Ts)
                f0[i] = 0;
            else if(T[i] > Tl)
                f0[i] = 1;
            else
                f0[i] = (T[i]-Ts)/(Tl-Ts);
        }
                
        //LIQUID MASS FRACTION - FUTURE TIME
        for(i=1; i<(n-1); i++)
        {
            if(T[i] < Ts)
                f[i] = f0[i];
            else if(T[i] > Tl)
                f[i] = f0[i];
            else
                f[i] = f0[i] + (1/(Tl-Ts))*(T[i]-Told[i]);
        }
                
        fprintf(sai, "%i        %e         %e \n", time, T[10], f0[10]);
        
        //EQUALS TEMPERATURES FIELD  
        for(i=0; i<n; i++)
            Told[i] = T[i];
        
        //COEFFICIENTS AND FORCE TERM
        for(i=1; i<(n-1); i++)
        {
            Kw = 0.5*(K[i-1]*K[i]); 
            Ke = 0.5*(K[i+1]*K[i]);
            Aw[i] = -Kw/(h*h);
            Ae[i] = -Ke/(h*h);
            Ap[i] = -(Ae[i] + Aw[i]) + ro*Cp[i]/dt;
            B[i] = -(ro*L/dt)*(f[i]-f0[i]) + ro*Cp0[i]*T[i]/dt;
        }
        
        //HETEROGENEOUS BOUNDARY CONDITION
        T[0] = 2*Tb - T[1];
        T[n-1] = 2*Tb - T[n-2];
                
        //PRINTING 
        if( (time == 0) || (fmod(time, limitTime/data) == 0) )
        {         
            FILE *fields;
            char title[15];
            sprintf(title, "FIELDS_%05i.dat", time);
            fields = fopen(title, "w");
            fprintf(fields, "variables = \"x\" , \"y\" , \"T\" , \"K\" , \"Cp0\" , \"f0\" \n");
            fprintf(fields, "zone i=%04d j=%04d f=point \n", nx, ny);
          
            vady = 0.0;
            for(j=1; j<(ny+ghost-1); j++)
            {
                vadx = 0.0;
                for(i=1; i<(nx+ghost-1); i++)
                {
                    fprintf(fields, "%e   %e   %e   %e   %e   %e \n", vadx, vady, T[i], K[i], Cp0[i], f0[i]);
                    vadx += h;
                }
                vady += h;
            }
            fclose(fields);
        }
        
        
        
        //SOLVER OF LINEAR SYSTEM    
        do
        {                      
            iteration++;
                      
            //S.O.R                       
            for(i=1; i<(n-1); i++)               
                T[i] = (omega/Ap[i])*(B[i] - (Ae[i]*T[i+1] + Aw[i]*T[i-1])) + (1.0 - omega)*T[i];
                    
            maximumResidue = 0.0;    
        
            //RESIDUE
            for(i=1; i<(n-1); i++)          
                R[i] = B[i] - (Ae[i]*T[i+1] + Aw[i]*T[i-1] + Ap[i]*T[i]);

            for(i=1; i<(n-1); i++)
            {
                if(fabs(R[i]) > maximumResidue)
                    maximumResidue = fabs(R[i]);
            }                    
            cout << "     MAXIMUM RESIDUE [" << iteration << "]: " << maximumResidue << endl;
        
            //cin >> pause;      
        }
        while(maximumResidue > eps);
        
        iteration = 0;
               
    }   
    
    fclose(sai);
    
    system("PAUSE");
    return EXIT_SUCCESS;
}
