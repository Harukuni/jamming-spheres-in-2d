/**********************************************************
O'Hern's algoirthm to generate the jammed configuration 
for gear model
*************************** *******************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "MT.h"
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;

const int N = 1024;
const double pi = M_PI;
const double err = 1.0E-16, err2 = 1.0E-15;
double L=1, dt = 0.05;

double X[N], Y[N], Vx[N], Vy[N];

#define MAX(a,b) (((a)>(b)) ? (a):(b))
#define MIN(a,b) (((a)<(b)) ? (a):(b))

/***************************************
Generatars of Random number 
*************************************/
int seed = 1;

#define  Random (double)rand()/RAND_MAX
int randi(int min, int max) {
  return genrand_int32()%(max-min+1)+min;
}

double uni(void) {
  return genrand_real3();
}

/*************************************************
For the particle diameters
*************************************************/

double S[N];

//sig2 > sig1 otherwise NL does not work.

double sig1 = 0.5, sig2 = 0.7;
double Sig(double x) {
  return sig1 + x*(sig2-sig1);
}

/***************************************
Interaction potential (harmonic)
**************************************/

double v(double h) {
  if(h<0) return 0.5*h*h;
  else return 0;
}
double dv(double h) {
  if(h<0) return h;
  else return 0;
}


/**********************************
For the boundary condition
*********************************/

void periodic(double &dx, double &dy) {
  while(dy > 0.5*L)  dy -= L;
  while(dy < -0.5*L) dy += L;
  while(dx > 0.5*L) dx -= L;
  while(dx < -0.5*L) dx += L;
  return;
}

double dis(int i, int j) {
  double dx = X[i]-X[j];
  double dy = Y[i]-Y[j];
  periodic(dx,dy);
  return sqrt(dx*dx + dy*dy);
}

void back_to_box(double &x, double &y) {
  while(x < 0) x += L;
  while(x > L) x-= L;

  while(y < 0) y += L;
  while(y > L) y -= L;

  return;
}

/***********************
 Gap function
**********************/

// h(i,j) = h(j,i) ?
double calc_h(int i, int j) {
  double dx = X[i]-X[j], dy = Y[i]-Y[j];
  periodic(dx,dy);
  double dr = sqrt(dx*dx + dy*dy);

  // gap function
  double s1 = S[i], s2 = S[j];

  return dr - s1 - s2;
}

// /*************************************************
// For the neighbor list
// *************************************************/


// //const int list_max = N;
// double r_NL = 0.1;
// double X0[N], Y0[N];
// int NL[N][N], count_NL[N];
// const int M_max = 1000;
// int Mx = N, My=N;
// int HEAD[M_max][M_max], LIST[N];
// double lx= L/Mx, ly=L/My;

// void make_list() {
//   // # of the cells
//   Mx = L/(2*sig2+r_NL);
//   My = L/(2*sig2+r_NL);
//   // the size of a cell
//   lx = L/Mx;
//   ly = L/My;
//   if(lx==0 or ly==0) {
//     printf("too small\n");
//     exit(-1);
//   }
//   // reset the HEADS
//   for(int i=0;i<Mx;i++) {
//     for(int j=0;j<My;j++) {
//       HEAD[i][j] = -1;
//     }
//   }
//   // construct the LISTs
//   for(int i=0;i<N;i++) {
//     int mx = X[i]/lx, my = Y[i]/ly;
//     // link to the previous occupant
//     LIST[i] = HEAD[mx][my];
//     // the last one goes to the header
//     HEAD[mx][my] = i;
//   }
//   return;
// }

// int bcx(int n) {
//   if(n==-1) return Mx-1;
//   else if(n==Mx) return 0;
//   else return n;
// }

// int bcy(int n) {
//   if(n==-1) return My-1;
//   else if(n==My) return 0;
//   else return n;
// }

// void update_NL(void) {
//   make_list();  
//   for(int i = 0; i < N; i++) {
//     X0[i] = X[i];
//     Y0[i] = Y[i];
//     count_NL[i] = 0;
//   }
//   for(int i=0;i<N;i++) {
//     //calcualte the cell number (mx,my,mz) 
//     int mx = X[i]/lx, my = Y[i]/ly;
//     //scan the nearest cells
//     for(int nx=mx-1;nx<mx+2;nx++)
//       for(int ny=my-1;ny<my+2;ny++) {
// 	int j = HEAD[bcx(nx)][bcy(ny)];
// 	//scan the particles in the cell (nx,ny,nz)
// 	while(j !=-1) {	  
// 	  if(i!=j and dis(i,j) < S[i]+sig2+r_NL) {
// 	    NL[i][count_NL[i]] = j;
// 	    count_NL[i]++;
// 	  }
// 	  j = LIST[j];
// 	}
//       }
//   }
//   return;
// }

// void check_NL() {
//   double dr_max = 0;
//   for(int i=0;i<N;i++) {

//     // relative distance
//     double dx = X[i]-X0[i], dy = Y[i]-Y0[i]; 
//     periodic(dx,dy);
//     double dr = sqrt(dx*dx + dy*dy);
//     if(dr > dr_max ) dr_max = dr;
//   }
//   if(2*dr_max > r_NL) {
//     update_NL();
//   }
//   return;
// }

const int list_max = N;
double r_NL = 0.1;
int NL[N][N];
int count_NL[N];


double X0[N], Y0[N];
void update_NL() {
  for(int i=0;i<N;i++) {
    for(int j=0;j<N;j++) {
      NL[i][j] = -1;
    }
  }

  for(int i=0;i<N;i++) {
    count_NL[i] = 0;
    for(int j=0;j<N;j++) {
      if(i==j) {continue;}

      // relative distance
      double dx = X[i]-X[j], dy = Y[i]-Y[j];
      periodic(dx,dy);
      double dr = sqrt(dx*dx + dy*dy);

      // add j-th particle to NL of i-th particle
      if(dr < S[i] + sig2 + r_NL) {
	NL[i][count_NL[i]] = j;
	count_NL[i]++;
      }
    }
  }

  // Meorize the positions when the NL was constructed
  for(int i=0;i<N;i++) {
    X0[i] = X[i];
    Y0[i] = Y[i];
  }

  return;
}

void check_NL() {
  double dr_max = 0;
  for(int i=0;i<N;i++) {
    // relative distance
    double dx = X[i]-X0[i], dy = Y[i]-Y0[i];
    periodic(dx,dy);
    double dr = sqrt(dx*dx + dy*dy);
    if(dr > dr_max ) dr_max = dr;
  }
  if(2*dr_max > r_NL) {
    update_NL();
  }
  return;
}



/**********************************
For the calculation of the forces
********************************/

double Fx[N], Fy[N];


double ene() {
  double sum = 0;
  for(int i=0;i<N;i++) {
    for(int n=0;n < count_NL[i];n++) {
      int j = NL[i][n];
      double h = calc_h(i,j);
      if(h>0) continue;
      sum += v(h);
    }
  }
  return 0.5*sum/N;
}



/**********************************************************
FIRE
**********************************************************/

double Fx1[N], Fy1[N];
double Fx2[N], Fy2[N];


//modified Euler
void MD2() {
  //calculate F(t)
  for(int i=0;i<N;i++) {
    Fx1[i] = 0;
    Fy1[i] = 0;

    for(int n=0;n< count_NL[i];n++) {
      int j = NL[i][n];

      double dx = X[i]-X[j], dy = Y[i]-Y[j];
      periodic(dx,dy);
      double dr = sqrt(dx*dx + dy*dy);

      // gap function
      double s1 = S[i], s2 = S[j];
      double h = dr - s1- s2;
      if(h>0) continue;

      // derivatives
      double r_x = dx/dr, r_y = dy/dr;
      double d_v = dv(h);

      // Force calculation
      Fx1[i] += -d_v*(r_x);
      Fy1[i] += -d_v*(r_y);
    }
  }

  //calculate x(t+dt)
  for(int i=0;i<N;i++) {
    X[i] += dt*Vx[i] + 0.5*dt*dt*Fx1[i];
    Y[i] += dt*Vy[i] + 0.5*dt*dt*Fy1[i];
    back_to_box(X[i],Y[i]);
  }
  check_NL();


  //calculate F(t+dt)
  for(int i=0;i<N;i++) {
    Fx2[i] = 0;
    Fy2[i] = 0;

    for(int n=0; n< count_NL[i];n++) {
      int j = NL[i][n];

      double dx = X[i]-X[j], dy = Y[i]-Y[j];
      periodic(dx,dy);
      double dr = sqrt(dx*dx + dy*dy);

      // gap function
      double s1 = S[i], s2 = S[j];
      double h = dr - s1 - s2;
      if(h>0) continue;

      // derivatives
      double r_x = dx/dr, r_y = dy/dr;
      double d_v = dv(h);

      // Force calculation
      Fx2[i] += -d_v*(r_x);
      Fy2[i] += -d_v*(r_y);

    }
  }

  //calculate forces and verocities
  for(int i=0;i<N;i++) {
    Fx[i] = 0.5*(Fx1[i] + Fx2[i]);
    Fy[i] = 0.5*(Fy1[i] + Fy2[i]);
    Vx[i] += dt*Fx[i];
    Vy[i] += dt*Fy[i];
  }
  return;
}


int Nmin = 5, Ns = 0;
double f_inc = 1.1, f_dec = 0.5, alpha_s = 0.1,
  f_alpha = 0.99, dt_max=0.1, alpha = 0.1;

void FIRE() {

  //MD2();
  MD2();
  double P = 0;
  for(int i=0;i<N;i++) {
    P += Fx[i]*Vx[i] + Fy[i]*Vy[i];
  }

  //calculate norm
  double Nf = 0, Nv = 0;
  for(int i=0;i<N;i++) {
    Nf += Fx[i]*Fx[i] + Fy[i]*Fy[i];
    Nv += Vx[i]*Vx[i] + Vy[i]*Vy[i];
  }
  Nf = sqrt(Nf), Nv = sqrt(Nv);

  //determine the direction
  for(int i=0;i<N;i++) {

    Vx[i] = (1-alpha)*Vx[i] + alpha*(Fx[i]/Nf)*Nv;
    Vy[i] = (1-alpha)*Vy[i] + alpha*(Fy[i]/Nf)*Nv;
  }

  if(P>0) {
    Ns++;
    if(Ns > Nmin) {
      dt = MIN(dt*f_inc, dt_max);
      alpha = alpha*f_alpha;
    }
  } else {
    dt = dt*f_dec;
    alpha = alpha_s;
    Ns = 0;
    for(int i=0;i<N;i++) {
      Vx[i] = 0;
      Vy[i] = 0;
    }
  }
  return;
}

int count_z() {
  int sum = 0;
  for(int i=0;i<N;i++) {
    for(int n=0;n< count_NL[i];n++) {
      int j = NL[i][n];
      if(calc_h(i,j) <0) sum ++;
    }
  }
  return sum;
}


double K() {
  double sum = 0;
  for(int i=0;i<N;i++) {
    sum += Fx[i]*Fx[i] + Fy[i]*Fy[i];
  }
  return sum/N;
}


int minimize() {
  //Reset parameters for FIRE
  dt_max = 0.2*sig1, dt = 0.01*sig1;
  Ns = 0, alpha = 0.1;

  int i = 0;
  int it = 11;
  while(1) {
    FIRE();
    if(K() < 1.0E-25 or ene() < 1.0E-16) break;
    // if(i>pow(2,it)) {
    //   double phi = 0.5*N*M_PI*(sig1*sig1 + sig2*sig2)/(L*L);
    //   printf("%d %f %e %e\n",i,phi,K(),ene());
    //   it++;
    // }
    i++;
  }
  return i;
}

/****************************
Calculation of contact number
*****************************/

bool Rattler[N];

int count_contact(int i ) {
  int sum = 0;
  if(Rattler[i]) return 0;

  for(int n=0;n< count_NL[i];n++) {
    int j = NL[i][n];
    if(Rattler[j]) continue;

    if(calc_h(i,j) <0) sum ++;
  }
  return sum;
}

int count_contact() {
  int sum = 0;
  for(int i=0;i<N;i++) {
    sum += count_contact(i);
  }
  return sum;
}

int count_rattlers() {
  int sum = 0;
  for(int i=0;i<N;i++) {
    if(Rattler[i]) sum++;
  }
  return sum;
}

void remove_rattlers() {

  for(int i=0;i<N;i++) Rattler[i] = false;

  while(1) {
    int z = count_contact();
    for(int i=0;i<N;i++) {
      if(Rattler[i])continue;
      double sum = 0;
      for(int n=0;n< count_NL[i];n++) {
	int j = NL[i][n];
	if(Rattler[j])continue;

	if(calc_h(i,j)<0) sum+=1.0;
      }
      if(sum < 2) Rattler[i] = true;
      else Rattler[i] = false;
    }
    int z_new = count_contact();
    if(z==z_new) break;
  }
  return;
}


/*********************************
Change the particl size for phi
******************************/

void set_phi(double phi) {
  // Set box size
  double L_old = L;

  double sum = 0;
  for(int i=0;i<N;i++) {
    sum += M_PI*S[i]*S[i];
  }

  L = sqrt(sum/phi);

  for(int i=0;i<N;i++) {
    X[i] = X[i]*L/L_old;
    Y[i] = Y[i]*L/L_old;
  }
  // Update neighbor list
  //check_NL();
  update_NL();
  return;
}


/*******************
Initial condition
*******************/

void init(double phi) {

  // Set box size
  double sum = 0;
  for(int i=0;i<N;i++) {
    double x = double(i)/(N-1);
    S[i] = Sig(x);
    sum += M_PI*S[i]*S[i];
  }
  L = sqrt(sum/phi);

  // Random initial condition
  for(int i=0;i<N;i++) {
    X[i] = L*uni();
    Y[i] = L*uni();
    Vx[i] = 0;
    Vy[i] = 0;
  }

  // Update neighbor list
  update_NL();
  return;
}


/*****************************
Find Jamming point
****************************/

const double dphi_end = 1.0E-15, e_end = 1.0E-16;

double phi_j2() {
  double phi = 0.1, dphi = 0.001;
  init(phi);
  int it = minimize();
  double e_old = ene();

  while(true) {
    phi += dphi;
    set_phi(phi);
    int it = minimize();
    double e_new = ene();

    if(e_new > e_end and e_new < 2*e_end)break;

    if(
       ((e_old < e_end) and (e_new > e_end))
              or
       ((e_old > e_end) and (e_new < e_end))
       ){dphi = -0.5*dphi;}

    e_old = e_new;

  }
  return phi;
}

/*****************************
Ajust pressure
****************************/

double pressure() {
  double sum = 0;
  for(int i=0;i<N;i++) {
    for(int n=0;n< count_NL[i];n++) {
      int j = NL[i][n];

      double dx = X[i]-X[j], dy = Y[i]-Y[j];
      periodic(dx,dy);
      double dr = sqrt(dx*dx + dy*dy);


      // gap function
      double s1 = S[i], s2 = S[j];

      double h = dr - s1- s2;

      //pressure
      if(h<0) {
	sum += -dr*dv(h);
      }
    }
  }
  return 0.5*sum/(L*L)/2;
}


double set_pressure(double p_end) {
  double dphi = p_end;
  //current pressure
  double p_old = pressure();
  
  //calculate phi
  double sum = 0;
  for(int i=0;i<N;i++) sum += pi*S[i]*S[i];
  double phi = sum/(L*L);
  
  while(true) {
    phi += dphi;
    set_phi(phi);
    int it = minimize();
    double p_new = pressure();

    if(p_new > p_end and p_new < 1.01*p_end)break;

    if(
       ((p_old < p_end) and (p_new > p_end))
              or
       ((p_old > p_end) and (p_new < p_end))
       ){dphi = -0.5*dphi;}

    //printf("phi=%f dphi=%e %e %e\n",phi,dphi,p_old,p_new);

    p_old = p_new;

  }
  return phi;
}


/***********************
Output
**********************/

void output(double p) {
  char filename[200];
  ofstream file;
  file<<setprecision(20);
  sprintf(filename,"N%dr%dp%e.txt",N,seed,p);
  file.open(filename);
  file << N << " " << seed << " " <<  L << " " << L << endl;
  for(int i=0;i<N;i++) {
    file << X[i] << " " << Y[i] << " " << S[i] << " " << Rattler[i] << endl;
  }
  file.close();
  return;
}

int main(int argc, char *argv[]) {
  seed = atoi(argv[1]);
  init_genrand(seed);
  double pmax = atof(argv[2]);
  

  // create jamming configuration
  double phij = phi_j2();
  double p = 0;
  remove_rattlers();
  int Nc = count_contact();
  int Nr = count_rattlers();
  double z = Nr<N ? Nc/double(N-Nr) : 0;
  printf("%e %f %f %e %d %d\n",p, phij, z, ene(), Nc, Nr);
  output(p);
  
  for(int i=0;i<6;i++) {
    for(int n=0;n<3;n++) {
      double p = pow(2,n)*pow(0.1,7-i);
      if(p > pmax)break;
      double phi = set_pressure(p);
      remove_rattlers();
      int Nc = count_contact();
      int Nr = count_rattlers();
      double z = Nr<N ? Nc/double(N-Nr) : 0;
      printf("%e %f %f %e %d %d\n",p, phi, z, ene(), Nc, Nr);
      output(p);
    }
  }
  return 0;
}
