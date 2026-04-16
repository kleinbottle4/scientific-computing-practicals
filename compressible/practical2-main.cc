// -*- origami-fold-style: triple-braces; -*-
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

struct solution {
  std::vector<double> u_old;
  std::vector<double> u_new;
  int npts;
  int ng; // The number of ghost points on each side.
  // The first non-ghost entry is u[ng].
  std::string scheme;
  std::string update_method;
  std::string initial_curve;
  double Dx, a, C;

  void set_initial_data(std::string curve)
  {
    for (int i = 1; i <= npts; i++) {
      double x = (i-1)*Dx - 1.0;
      if (curve == "gaussian"){
	u_old[i] = exp(-8*x*x);
      } else if (curve == "tophat"){
	u_old[i] = (0.25 < x && x < 0.75 ?  2 : 1);
      } else {
	u_old[i] = 0;
      }
    }
    // Ghost points
    u_old[0] = u_old[1];
    u_old[npts + 1] = u_old[npts];
  }

  void print()
  {
    for (int i=1; i <= npts; i++) {
      std::cout << (i-1.0)*Dx - 1.0 << '\t' << u[i] << std::endl;
    }
  }

  void LF_update()
  {
    // update the solution to the next time step using
    // finite differences lax-friedrichs
    for (int i = 1; i <= npts; i++){
      u_new[i] = 0.5*(1.0-C)*u[i+1] + 0.5*(1.0+C)*u[i-1];
    }
    // Ghost points:
    u_new[ng - 1] = u[ng + 0];
    u_new[ng + npts] = u[ng + npts - 1];

    //
    u = u_new;
  }

  void LW_update()
  {
    // Update the solution to the next time-step using the
    // Lax-Wendroff method.
    for (int i=1; i <= npts; i++) {
      u_new[i] = (1-C*C)*u_old[i] + 0.5*C*(1+C)*u_old[i-1] - 0.5*C*(1-C)*u_old[i+1];
    }

    // Ghost points:
    u_new[0] = u_new[1];
    u_new[u_new.size()] = u_new[u_new.size() - 1];
  }

  void warming_beam_update()
  {

  }
};


int main()
{
  solution u; 

  // get parameters from user
  u.npts = 200;
  u.Dx = 0.01;
  u.a = 1;
  u.C = 0.8;
  double Dt = C*Dx / a;
  double T = 1.0;

  u.set_initial_data("gaussian");

  u.print();

  for (double t = 0.0; t < T; t += Dt){
    u.update();
  }

  // print the data
  u.print();

  return 0;
}
