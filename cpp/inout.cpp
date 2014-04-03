// usage: ./inout nx ny Cx Cy nt it 1>...:in 2>...:out

#include "listings.hpp"
#include <iostream>

using namespace std;
using namespace blitz;

ostream& operator<<(ostream& os, const Array<real_t, 2>& arr)
{
  for (int i=arr.lbound(0); i<=arr.ubound(0); ++i) 
  {
    os << arr(i, 0);
    for (int j=arr.lbound(1)+1; j<=arr.ubound(1); ++j) 
    {
      os << "\t" << arr(i, j);
    }
    os << "\n";
  }
  return os;
}

template <int it> 
void inout(char **argv)
{
  int 
    nx = atoi(argv[1]), 
    ny = atoi(argv[2]),
    nt = atoi(argv[5]);
  real_t
    Cx = atof(argv[3]),
    Cy = atof(argv[4]);

  solver_mpdata<it, cyclic<0>, cyclic<1>> slv(nx, ny); 
  {
    firstIndex i;
    secondIndex j;
    slv.state() = exp(
      -sqr(i-nx/2.) / (2.*pow(nx/10., 2)) 
      -sqr(j-ny/2.) / (2.*pow(ny/10., 2)) 
    );
  }
  slv.courant(0) = Cx; 
  slv.courant(1) = Cy; 

  cout << slv.state();
  slv.solve(nt);
  cerr << slv.state();
}

int main(int argc, char **argv)
{
  if (argc != (6 + 1))
  {
    cerr << "expecting 6 arguments (nx, ny, Cx, Cy, nt, it)" << endl;
    throw; 
  }

  int it = atoi(argv[6]);

  try
  {
    switch (it)
    {
      case 1: inout<1>(argv); break;
      case 2: inout<2>(argv); break;
      case 3: inout<3>(argv); break;
      case 4: inout<4>(argv); break;
      case 5: inout<5>(argv); break;
      case 6: inout<6>(argv); break;
      default: throw;
    }
  }
  catch (...)
  {
    exit(1);
  }
}
