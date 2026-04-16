// file practical3/main.cc: To solve Burger's equation with the Finite
// Volume Method under various schemes.

// Compile with -D_GLIBCXX_DEBUG to get bounds checking for vectors

#include <iostream>
#include <vector>

struct solution {
  /*
    The mid-point of each cell is (i-0.5)*Dx.
  */
  double Dx;
  std::vector<double> u;
  std::vector<double> u_new;
  std::vector<double> f; // store f_{i-1/2} in f[i]
  int ncells;
  double Dt;

  std::vector<std::vector<double>> save_data;

  solution(int number_of_cells)// {{{
  {
    ncells = number_of_cells;
    u.resize(ncells+2);
    u_new.resize(ncells+2);
    f.resize(ncells+1);
  }
  // }}}

  void set_initial_data()// {{{
  {
    for (int i = 0; i < ncells; i++){
      double x = (i+0.5)*Dx;
      u[1 + i] = (x < 0.5 ? 2 : 1);
    }
    save_data.push_back(u);

    // Set initial ghost points
    u[0] = u[1];
    u[ncells + 1] = u[ncells];
  }
  // }}}

  void set_initial_flux()// {{{
  {
    // Set f_{-1/2} to 
    f[0] = 0.5*u[0]*u[0];

    f[ncells] = 0.5*u[1 + ncells]*u[1 + ncells];
  }
  // }}}
  double flux(double u) const
  {
    // The true flux, for Burger's equation.
    return 0.5*u*u;
  }


  double LF_flux(int i) const
  {
    // f_{i+1/2}^{LF} // Uses u_i^n and u_{i+1}^n
    return 0.5*(Dx/Dt)*(u[1+i] - u[i+2]) +
      0.5*(flux(u[i+1]) + flux(u[i+2]));
  }

  double RI_flux(int i) const
  {
    double z; // u_{i + 1/2}^{n + 1/2}
    z = 0.5*(u[i] + u[i+1]) - 0.5*(Dt/Dx)*(flux(u[i+1])-flux(u[i]));
    return flux(z);
  }

  double FORCE_flux(int i) const
  {
    return 0.5*(LF_flux(i) + RI_flux(i));
  }

  double numerical_flux(int i) const
  {
    return FORCE_flux(i);
  }

  void update()
  {
    for (int i = 0; i < ncells; i++){
      // f_{i+1/2} = 
      f[i+1] = numerical_flux(i+1);
      // u_{i}^{n+1} = u_{i}^{n} - (Dt/Dx)*(f_{i+1/2} - f_{i-1/2})
      u_new[1+i] = u[1+i] - (Dt/Dx)*(f[i+1]-f[i]);
    }
    u = u_new;
    // The ghost points
    u[0] = u[1];
    u[1+ncells] = u[ncells];
    save_data.push_back(u);
  }

  void print()
  {
    // Print the data in columns, with the first column being the x coordinates
    // Don't print the ghost cells.
    for (int i = 0; i < ncells; i++){
      double x = (i+0.5)*Dx;
      std::cout << x << '\t';
      for (size_t j = 0; j < save_data.size(); j++){
	std::cout << save_data[j][1+i] << '\t';
      }
      std::cout << std::endl;
    }
  }
};

int main()
{
  int ncells = 100;
  solution u(ncells);
  u.Dx = 0.01;
  double Dt = 0.001;
  u.Dt = Dt;
  double T = 0.2;

  u.set_initial_data();
  u.set_initial_flux();
		
  double t = 0;

  while (t < T){
    t += Dt;
    u.update();
  }

  u.print();

  return 0; 
}
