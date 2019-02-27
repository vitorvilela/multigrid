//FRONTEIRAS

#ifndef fronteirasH
#define fronteirasH


void Fronteira_P(double ***P, double h, int nx, int ny, int nz, char fron, int sit);

void Fronteira_U(double ***U, double h, int nx, int ny, int nz, char fron);

void Fronteira_V(double ***V, double h, int nx, int ny, int nz, char fron);

void Fronteira_W(double ***W, double h, int nx, int ny, int nz, char fron);

void Periodica_P(double ***P, double h, int nx, int ny, int nz);

void Periodica_U(double ***U, double h, int nx, int ny, int nz);

void Periodica_V(double ***V, double h, int nx, int ny, int nz);

void Periodica_W(double ***W, double h, int nx, int ny, int nz);

void Dirichlet(int nx, int ny, int nz, double ***U, double ***V, double ***W, double ***Uex, double ***Vex, double ***Wex);

void Update_vel(int nx, int ny, int nz, double ***Ue, double ***Ve, double ***We, double ***Uex, double ***Vex, double ***Wex);

void Cavidade_vel(int nx, int ny, int nz, double ***U, double ***V, double ***W, double Ut);

void Update_vel_cavidade(int nx, int ny, int nz, double ***U, double ***V, double ***W);

void Cavidade_pressao(int nx, int ny, int nz, double ***Pe, double h, int sit);


#endif
