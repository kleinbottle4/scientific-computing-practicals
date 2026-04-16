#include <iostream>
#include <map>
#include <cassert>

double flux(double u)
{
	// The true flux, for Burger's equation.
//	return 0.5*u*u;
	return u;
}

double numerical_flux(double i, double Dx, double Dt, std::map<int, double> &u)
{
	//Lux Friedrichs
	// f_{i+1/2}^{LF}
	// u is u_i^n
	// v is u_{i+1}^n
	return 0.5*(Dx/Dt)*(u[i] - u[i+1]) + 0.5*(flux(u[i]) + flux(u[i+1]));
}

int main()
{
	double Dx = 0.01;
	std::map<int, double> u;
	std::map<int, double> u_new;
	std::map<double, double> f;


	// Set initial data
	for (int i = 0; i <= 99; i++){
		double x = i/100.0 + Dx/2.0;
		u[i] = (x < 0.5 ? 2 : 1);
		std::cout << x << '\t' << u[i] << std::endl;
	}

	// Set initial ghost points
	// Transmissive boundary
	u[-1] = u[0];
	u[100] = u[99];

	std::cout << std::endl;

	double t = 0;

	int counter = 1;
	while (t < 0.2){
		double Dt = 0.0001;
		t += Dt;
		for (int i = 0; i <= 99; i++){
			double x = i/100.0 + Dx/2;
			// f_{i-1/2}
			f[i - 0.5] = numerical_flux(i, Dx, Dt, u);
			// f_{i+1/2} = 
			f[i + 0.5] = numerical_flux(i, Dx, Dt, u);
			// u_{i}^{n+1} = u_{i}^{n} - (Dt/Dx)*(f_{i+1/2} - f_{i-1/2})
			u_new[i] = u[i] - (Dt / Dx) * (f[i+0.5] - f[i-0.5]);
			if (counter % 10 == 0)
				std::cout << x << '\t' << u_new[i] << std::endl;
		}
		counter++;
		u = u_new;
		// The ghost points
		u[-1] = u[0];
		u[100] = u[99];
		std::cout << std::endl;
	}
	return 0;
}
