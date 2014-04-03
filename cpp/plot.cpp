//listing20
#include "listings.hpp"
#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

enum {x, y};

template <class T>
void setup(T &solver, int n[2]) 
{
  blitz::firstIndex i;
  blitz::secondIndex j;
  solver.state() = exp(
    -sqr((.5+i)-n[x]/2.) / (2*pow(n[x]/10., 2))
    -sqr((.5+j)-n[y]/2.) / (2*pow(n[y]/10., 2))
  );  
  solver.courant(x) = -.5; 
  solver.courant(y) = -.25;
}

template <class T>
void plot(T &solver, Gnuplot &gp)
{
  gp << "splot '-' binary" 
     << gp.binfmt(solver.state())
     << " origin=(.5,.5,-1)"
     << " with image notitle"
     << ", '-' binary" 
     << gp.binfmt(solver.state())
     << " origin=(.5,.5,0)"
     << " with lines notitle\n";
  gp.sendBinary(solver.state().copy());
  gp.sendBinary(solver.state().copy());
}

int main() 
{
  int n[] = {24, 24}, nt = 75;
  Gnuplot gp;
  gp << "set term pdf size 10cm, 30cm\n" 
     << "set output 'figure.pdf'\n"     
     << "set multiplot layout 4,1\n" 
     << "set border 4095\n"
     << "set xtics out\n"
     << "set ytics out\n"
     << "unset ztics\n"    
     << "set xlabel 'x/dx'\n"
     << "set ylabel 'y/dy'\n"
     << "set xrange [0:" << n[x] << "]\n"   
     << "set yrange [0:" << n[y] << "]\n"   
     << "set zrange [-1:1]\n"   
     << "set cbrange [-.025:1.025]\n"     
     << "set palette maxcolors 42\n";
  {
    solver_donorcell<cyclic<x>, cyclic<y>> 
      slv(n[x], n[y]);
    setup(slv, n);
    gp << "set title 't/dt=0'\n";
    plot(slv, gp);
    slv.solve(nt);
    gp << "set title 'donorcell t/dt="<<nt<<"'\n";
    plot(slv, gp);
  } 
  {
    const int it = 2;
    solver_mpdata<it, cyclic<x>, cyclic<y>> 
      slv(n[x], n[y]); 
    setup(slv, n); 
    slv.solve(nt);
    gp << "set title 'mpdata<" << it << "> "
       << "t/dt=" << nt << "'\n";
    plot(slv, gp);
  } 
  {
    const int it = 44;
    solver_mpdata<it, cyclic<x>, cyclic<y>> 
      slv(n[x], n[y]); 
    setup(slv, n); 
    slv.solve(nt); 
    gp << "set title 'mpdata<" << it << "> "
       << "t/dt=" << nt << "'\n";
    plot(slv, gp);
  }
}
//listing21
