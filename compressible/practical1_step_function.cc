#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

#define PI 3.14159265358979323846264

void print_results(std::ofstream & os, const std::vector<double> &u, double Dx,
		int num_pts, int num_ghost_pts)
{
	for (int i = 0; i < num_pts; i++){
		os << i*Dx << " " << u.at(i + num_ghost_pts) << ' ' << std::endl;
	}
}


int main()
{
	const int npts = 100;
	std::vector<double> old_u(npts + 2);
	std::vector<double> new_u(npts + 2);

	const double Dx = 0.01;
	const double Dt = 0.01;
	const double a = 1;
	const double gamma = Dx + a*Dt;

	for (int i = 1; i <= npts; i++){
		double x = double(i-1) * Dx;
//		old_u.at(i) = sin(2*PI*x);
		old_u.at(i) = ((0.25 < x && x < 0.75) ? 1 : 0);
	}
	// Ghost points (periodic boundary conditions)
	old_u.at(0) = old_u.at(npts);
	old_u.at(npts+1) = old_u.at(1);

	std::ofstream output_file("advection_data.dat", std::ios_base::app);
	print_results(output_file, old_u, Dx, npts, 1);

	double t = 0;
	while (t <= 1.0){
		t += Dt;
		for (int i = 1; i <= npts; i++)
			new_u.at(i) = (Dx*old_u.at(i) + a*Dt*new_u.at(i-1)) / gamma;
		// ghost points (periodic boundary conditions)
		new_u.at(0) = new_u.at(npts);
		new_u.at(npts + 1) = new_u.at(1);

		std::vector<double> tmp = old_u;
		old_u = new_u;
		new_u = tmp;
	}

	print_results(output_file, new_u, Dx, npts, 1);
	output_file.close();

	return 0;
}
