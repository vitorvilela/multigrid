//solvers.h

#ifndef solversH
#define solversH

void nested(double ****Pex, double ****E, double ****e, double ****R, double ****rhs, double ****z1, double ****z2, double ****R1, double ****R2, double ****p, double ****W, double ****Ap, double ****Ae, double ****Aw, double ****As, double ****An, double ****Au, double ****Ad, double ****Mp, double ****Me, double ****Mw, double ****Ms, double ****Mn, double ****Mu, double ****Md, double ****B, double ****P, int *nx, int *ny, int *nz, double *h, int niveis, int fun, char fron, double L, double D, double S, int *itr_u, int *itr_d, double eps, double *w);

int MG(int niveis, double ****e, double ****R, double ****rhs, double ***P, double ***B, double ****Ap, double ****Ae, double ****Aw, double ****An, double ****As, double ****Au, double ****Ad, int *itr_u, int *itr_d, int *nx, int *ny, int *nz, double D, double L, double S, double *h, double *w, double eps, char fron);

void MGgM(double ***Pex, double ***E, double ****p, double ****R, double ****K, double ****S, double ****Ap, double ****Ae, double ****Aw, double ****As, double ****An, double ****Au, double ****Ad, double ****B, double ****P, int *nx, int *ny, int *nz, double *h, int niveis, int fun, char fron, double L, double D, double s, double eps, char pre);

void sor(double ***Pe, double ***B, double ***Ap, double ***Ae, double ***Aw, double ***An, double ***As, double ***Au, double ***Ad, int nx, int ny, int nz, double D, double L, double S, double h, double w, int itr, int sit, char fron);

void cg(double ***Pex, double ***P, double ***B, double ***E, double ***z1, double ***z2, double ***R1, double ***R2, double ***p, double ***W, double ****R, double ****rhs, double ****e, double ****Ap, double ****Ae, double ****Aw, double ****An, double ****As, double ****Au, double ****Ad, double ****Mp, double ****Me, double ****Mw, double ****Mn, double ****Ms, double ****Mu, double ****Md, double ****Zero, int *nx, int *ny, int *nz, double D, double L, double S, double *h, double *w, int *itr_d, int *itr_u, double eps, int fun, char fron, char pre, int niveis);

void jacobi(double ***P, double ***B, double ***Ap, int nx, int ny, int nz);

void gradiente(double ***Pex, double ***E, double ***P, double ***B, double ***R, double ***W, double ***Ap, double ***Aw, double ***Ae, double ***An, double ***As, double ***Au, double ***Ad, int nx, int ny, int nz, double D, double L, double S, double h, double eps, int fun, char fron);

#endif

