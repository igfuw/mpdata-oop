#include "listings.hpp"
#include <fstream>
#include <ios>

using namespace std;

inline arr_t read_file(const string &file, int nx, int ny)
{
  arr_t array(nx, ny);
  stringstream io;
  io << "(0," << nx-1 << ") x (0," << ny-1 << ")" << endl << "[";
  io << ifstream(file, ios::binary).rdbuf();
  io << "]";
  io >> array;
  return array;
}

template <int it> 
void test(char **argv)
{
  int 
    nt = atoi(argv[5]),
    nx = atoi(argv[1]), 
    ny = atoi(argv[2]),
    dec = atoi(argv[9]);
  real_t
    Cx = atof(argv[3]),
    Cy = atof(argv[4]);
  string 
    fin(argv[7]),
    fout(argv[8]);

  solver_mpdata<it, cyclic<0>, cyclic<1>> slv(nx, ny); 
  slv.state() = read_file(fin, nx, ny);
  slv.courant(0) = Cx; 
  slv.courant(1) = Cy; 
  slv.solve(nt);
  if (max(abs(slv.state() - read_file(fout, nx, ny))) >= .5 * pow(10, -dec)) 
  {
    cerr << slv.state() << endl;
    throw;
  }
}

int main(int argc, char **argv)
{
  if (argc != (9 + 1))
  {
    cerr << "expecting 9 arguments (nx, ny, Cx, Cy, nt, it, f.in, f.out, dec)" << endl;
    throw; 
  }

  int it = atoi(argv[6]);

  try
  {
    switch (it)
    {
      case 1: test<1>(argv); break;
      case 2: test<2>(argv); break;
      case 3: test<3>(argv); break;
      case 4: test<4>(argv); break;
      case 5: test<5>(argv); break;
      case 6: test<6>(argv); break;
      default: throw;
    }
  }
  catch (...)
  {
    exit(1);
  }
}
