from listings import *
import sys

def read_file(fname, nx, ny):
  tmp = numpy.empty((nx, ny), real_t)
  with open(fname, 'r') as f:
    x = 0
    for line in f:
      tmp[x,:] = numpy.fromstring(line, dtype=real_t, sep='\t')
      x += 1
  assert(x == nx)
  return tmp

def main():
  if (len(sys.argv) != (9 + 1)) : 
    raise Exception('expecting 9 arguments (nx, ny, Cx, Cy, nt, it, f.in, f.out, dec)')

  nx = int(sys.argv[1])
  ny = int(sys.argv[2])
  Cx = float(sys.argv[3])
  Cy = float(sys.argv[4])
  nt = int(sys.argv[5])
  it = int(sys.argv[6])
  fname = sys.argv[7]
  fout = sys.argv[8]
  dec = int(sys.argv[9])

#listing19
  slv = solver_mpdata(it, cyclic, cyclic, nx, ny)
  slv.state()[:] = read_file(fname, nx, ny)
  slv.courant(0)[:] = Cx
  slv.courant(1)[:] = Cy
  slv.solve(nt)
#listing20

  diff = numpy.amax(abs(slv.state() - read_file(fout, nx, ny)))
  print "diff=", diff
  if (diff >= .5 * pow(10, -dec)): 
    print slv.state().dtype, slv.state()
    tmp = read_file(fout, nx, ny)
    print tmp.dtype, tmp
    print numpy.max(numpy.abs(slv.state() - tmp))
    raise Exception()

if __name__ == "__main__":
    main()
