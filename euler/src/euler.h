//EULER

#ifndef eulerH
#define eulerH


void ue_laminar(int nx, int ny, int nz, double dt, double ***Ro, double mi, double ***Ue, 
                double ***U, double ***V, double ***W, double ***P, double ***Fx, double h);
                
void ve_laminar(int nx, int ny, int nz, double dt, double ***Ro, double mi, double ***Ve, 
                double ***U, double ***V, double ***W, double ***P, double ***Fy, double h);
                
void we_laminar(int nx, int ny, int nz, double dt, double ***Ro, double mi, double ***We, 
                double ***U, double ***V, double ***W, double ***P, double ***Fz, double h);
                
void f(int nx, int ny, int nz, double ***B, double ***Ue, double ***Ve, double ***We, double dt, double h);

void p(int nx, int ny, int nz, double ***P, double ***Pe);

void u(int nx, int ny, int nz, double dt, double ***Ro, double ***U, double ***Ue, double ***Pe, double h);

void v(int nx, int ny, int nz, double dt, double ***Ro, double ***V, double ***Ve, double ***Pe, double h);

void w(int nx, int ny, int nz, double dt, double ***Ro, double ***W, double ***We, double ***Pe, double h);
    
void continuidade(double ***U, double ***V, double ***W, int nx, int ny, int nz, double h, int tempo);


#endif


