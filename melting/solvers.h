//SOLVERS


#ifndef solversH
#define solversH

void redblack(double **P, double **B, double **Ap, double **Ae, double **Aw, double **An, double **As, int nx, int ny, double L, double D, double h, double w, int itr, int sit, int fun, char fron);

void sor(double **T, double **K, double **B, double **Ap, double **Ae, double **Aw, double **An, double **As, int nx, int ny, double L, double D, double h, double w, int itr);
void jacobi(double **P, double **B, double **Ap, int nx, int ny);

void cg(double **R1, double **R2, double **z1, double **z2, double **p, double **W, double **P, double **B, double ***Ap, double ***Ae, double ***Aw, double ***As, double ***An, int *nx, int *ny, double L, double D, double *h, int fun, char fron, int niveis, double eps, char pre, int *itr_u, int *itr_d, double *w, double ***e, double ***rhs, double ***R, double ***Zero);

void multigrid(double ***e, double ***R, double ***rhs, double **P, double **B, double ***Ap, double ***Ae, double ***Aw, double ***An, double ***As, int *itr_u, int *itr_d, int *nx, int *ny, double L, double D, double *h, double *w, int fun, char fron, double eps, int niveis);

void nested(double ***e, double ***R, double ***rhs, double ***Ap, double ***Ae, double ***Aw, double ***As, double ***An, double ***Pex, double ***B, double ***P, int *nx, int *ny, double *h, int niveis, int fun, char fron, double L, double D, int *itr_u, int *itr_d, double eps, double *w);

void MGgM(double ***p, double ***R, double ***K, double ***S, double ***Ap, double ***Ae, double ***Aw, double ***As, double ***An, double ***B, double ***P, int *nx, int *ny, double *h, int niveis, int fun, char fron, double L, double D, double eps, char pre);








#endif
