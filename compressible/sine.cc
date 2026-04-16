#include <iostream>
#include <map>
#include <cmath>
#include <fstream>
#include <cassert>

#define PI 3.14159265358979323846264


/* partial u / partial t + a partial u / partial x = 0 */

int main(int argc, const char **argv)
{
  if (argc < 5){
    std::cout << "Usage. a.out <a> <scheme> <T> <Delta T>.\n";
    return 1;
  }

  std::map<int, double> u;
  std::map<int, double> new_u;

  const double Dx = 0.01;

  const double a = atof(argv[1]);
  const int scheme = atoi(argv[2]);
  const double T = atof(argv[3]);
  const double Dt = atof(argv[4]);

  std::cout << "Received a = " << a << ".\n";
  std::cout << "Received scheme = " << scheme << ".\n";
  std::cout << "Received T = " << T << ".\n";
  std::cout << "Received Delta t = " << Dt << ".\n";

  assert(scheme == 0 || scheme == 1 || scheme == -1);


  for (int i = 0; i <= 100; i++){
    double x = double(i) * Dx;
    u[i] = sin(2*PI*x);
  }
  // Ghost points (periodic boundary conditions)
  u[-1] = u[100];
  u[101] = u[0];

  std::ofstream output_file("advection_data_sine.dat");

  for (int i = 0; i <= 100; i++)
    output_file << i*Dx << '\t' << u[i] << std::endl;

  double t = 0;
  while (t < T) {
    t += Dt;
    for (int i = 0; i <= 100; i++){
      switch(scheme){
      case 1:
	// forward spatial derivative
	new_u[i] = u[i] - a*(Dt/Dx)*u[i+1] + a*(Dt/Dx)*u[i]; break;
      case -1:
	// backwards spatial derivative
	new_u[i] = u[i] - a*(Dt/Dx)*u[i] + a*(Dt/Dx)*u[i-1];
	break;
      case 0:
	// centered spatial derivative
	new_u[i] = u[i] - 0.5*a*(Dt/Dx)*u[i+1] + 0.5*a*(Dt/Dx)*u[i-1];
	break;
      }
    }
    // ghost points (periodic boundary conditions)
    new_u[-1] = new_u[100];
    new_u[101] = new_u[0];

    u = new_u;
  }

  for (int i = 0; i <= 100; i++)
    output_file << i*Dx << '\t' << new_u[i] << std::endl;

  output_file.close();

  return 0;
}
