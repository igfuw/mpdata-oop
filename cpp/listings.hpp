//listing00
// code licensed under the terms of GNU GPL v3
// copyright holder: University of Warsaw
//listing01
using real_t = double;
//listing02
#include <blitz/array.h>
using arr_t = blitz::Array<real_t, 2>;
using rng_t = blitz::Range;
using idx_t = blitz::RectDomain<2>;
//listing03
#define return_macro(expr) \
  -> decltype(safeToReturn(expr)) \
{ return safeToReturn(expr); } 
//listing04
#include <boost/ptr_container/ptr_vector.hpp>
struct arrvec_t : boost::ptr_vector<arr_t> 
{
  const arr_t &operator[](const int i) const 
  {   
    return this->at((i + this->size()) % this->size()); 
  }
};
//listing05
struct hlf_t {} h;

inline rng_t operator+(const rng_t &i, const hlf_t &) 
{ 
  return i; 
} 

inline rng_t operator-(const rng_t &i, const hlf_t &) 
{ 
  return i-1; 
}
//listing06
template<int d> 
inline idx_t pi(const rng_t &i, const rng_t &j);

template<>
inline idx_t pi<0>(const rng_t &i, const rng_t &j) 
{
  return idx_t({i,j});
};

template<>
inline idx_t pi<1>(const rng_t &j, const rng_t &i) 
{
  return idx_t({i,j});
}; 
//listing07
template<class T1, class T2, class T3> 
inline auto F(
  const T1 &psi_l, const T2 &psi_r, const T3 &C
) return_macro(
  (
    (C + abs(C)) * psi_l + 
    (C - abs(C)) * psi_r
  ) / 2
)
//listing08
template<int d>  
inline auto donorcell( 
  const arr_t &psi, const arr_t &C, 
  const rng_t &i, const rng_t &j
) return_macro(
  F(
    psi(pi<d>(i,   j)), 
    psi(pi<d>(i+1, j)), 
      C(pi<d>(i+h, j))
  ) -
  F(
    psi(pi<d>(i-1, j)), 
    psi(pi<d>(i,   j)), 
      C(pi<d>(i-h, j))
  )
)
//listing09
void donorcell_op(
  const arrvec_t &psi, const int n,
  const arrvec_t &C, 
  const rng_t &i, const rng_t &j
) { 
  psi[n+1](i,j) = psi[n](i,j) - (
    donorcell<0>(psi[n], C[0], i, j) +
    donorcell<1>(psi[n], C[1], j, i)
  ); 
}
//listing10
template<class nom_t, class den_t>
inline auto mpdata_frac(
  const nom_t &nom, const den_t &den
) return_macro(
  where(den > 0, nom / den, 0)
) 
//listing11
template<int d>
inline auto mpdata_A(const arr_t &psi, 
  const rng_t &i, const rng_t &j
) return_macro(
  mpdata_frac(
    psi(pi<d>(i+1, j)) - psi(pi<d>(i,j)),
    psi(pi<d>(i+1, j)) + psi(pi<d>(i,j))
  ) 
) 
//listing12
template<int d>
inline auto mpdata_B(const arr_t &psi, 
  const rng_t &i, const rng_t &j
) return_macro(
 mpdata_frac(
    psi(pi<d>(i+1, j+1)) + psi(pi<d>(i, j+1)) -
    psi(pi<d>(i+1, j-1)) - psi(pi<d>(i, j-1)),
    psi(pi<d>(i+1, j+1)) + psi(pi<d>(i, j+1)) +
    psi(pi<d>(i+1, j-1)) + psi(pi<d>(i, j-1))
  ) / 2
)
//listing13
template<int d>
inline auto mpdata_C_bar(
  const arr_t &C, 
  const rng_t &i, 
  const rng_t &j
) return_macro(
  (
    C(pi<d>(i+1, j+h)) + C(pi<d>(i, j+h)) +
    C(pi<d>(i+1, j-h)) + C(pi<d>(i, j-h)) 
  ) / 4
)
//listing14
template<int d>
inline auto mpdata_C_adf(
  const arr_t &psi, 
  const rng_t &i, const rng_t &j,
  const arrvec_t &C
) return_macro(
  abs(C[d](pi<d>(i+h, j))) 
  * (1 - abs(C[d](pi<d>(i+h, j)))) 
  * mpdata_A<d>(psi, i, j) 
  - C[d](pi<d>(i+h, j)) 
  * mpdata_C_bar<d>(C[d-1], i, j)
  * mpdata_B<d>(psi, i, j)
) 
//listing15
template<class n_t>
inline rng_t ext(const rng_t &r, const n_t &n) { 
  return rng_t(
    (r - n).first(), 
    (r + n).last()
  ); 
} 
//listing16
template<class bcx_t, class bcy_t>
struct solver
{
  // member fields
  arrvec_t psi, C;
  int n, hlo;
  rng_t i, j;
  bcx_t bcx;
  bcy_t bcy;

  // ctor
  solver(int nx, int ny, int hlo) :
    hlo(hlo),
    n(0), 
    i(0, nx-1), 
    j(0, ny-1),  
    bcx(i, j, hlo), 
    bcy(j, i, hlo)
  {
    for (int l = 0; l < 2; ++l) 
      psi.push_back(new arr_t(ext(i, hlo), ext(j, hlo)));
    C.push_back(new arr_t(ext(i, h), ext(j, hlo)));
    C.push_back(new arr_t(ext(i, hlo), ext(j, h)));
  }

  // accessor methods
  arr_t state() {
    return psi[n](i,j).reindex({0,0});
  }

  arr_t courant(int d) 
  { 
    return C[d]; 
  }

  // helper methods invoked by solve()
  virtual void advop() = 0;

  void cycle() 
  { 
    n = (n + 1) % 2 - 2; 
  }

  // integration logic
  void solve(const int nt) 
  {
    for (int t = 0; t < nt; ++t) 
    {
      bcx.fill_halos(psi[n], ext(j, hlo));
      bcy.fill_halos(psi[n], ext(i, hlo));
      advop();
      cycle();
    }
  }
};
//listing17
template<int d>
struct cyclic
{
  // member fields
  rng_t left_halo, rght_halo;
  rng_t left_edge, rght_edge;;

  // ctor
  cyclic(
    const rng_t &i, const rng_t &j, int hlo
  ) :
    left_halo(i.first()-hlo, i.first()-1),
    rght_edge(i.last()-hlo+1, i.last()  ),
    rght_halo(i.last()+1, i.last()+hlo  ),
    left_edge(i.first(), i.first()+hlo-1)
  {} 

  // method invoked by the solver
  void fill_halos(const arr_t &a, const rng_t &j)
  {
    a(pi<d>(left_halo, j)) = a(pi<d>(rght_edge, j));     
    a(pi<d>(rght_halo, j)) = a(pi<d>(left_edge, j));     
  }
};
//listing18
template<class bcx_t, class bcy_t>
struct solver_donorcell : solver<bcx_t, bcy_t> 
{
  solver_donorcell(int nx, int ny) :
    solver<bcx_t, bcy_t>(nx, ny, 1)
  {}  

  void advop()
  {
    donorcell_op(
      this->psi, this->n, this->C, 
      this->i, this->j
    );
  }
};
//listing19
template<int n_iters, class bcx_t, class bcy_t>
struct solver_mpdata : solver<bcx_t, bcy_t>
{
  // member fields
  static const int n_tmp = n_iters > 2 ? 2 : 1;
  arrvec_t tmp[n_tmp];
  rng_t im, jm;

  // ctor
  solver_mpdata(int nx, int ny) : 
    solver<bcx_t, bcy_t>(nx, ny, 1), 
    im(this->i.first() - 1, this->i.last()),
    jm(this->j.first() - 1, this->j.last())
  {
    for (int n = 0; n < n_tmp; ++n)
    {
      tmp[n].push_back(new arr_t(
        this->C[0].domain()[0], this->C[0].domain()[1])
      );
      tmp[n].push_back(new arr_t(
        this->C[1].domain()[0], this->C[1].domain()[1])
      );
    }
  }

  // method invoked by the solver
  void advop()
  {
    for (int step = 0; step < n_iters; ++step) 
    {
      if (step == 0) 
        donorcell_op(
          this->psi, this->n, this->C, this->i, this->j
        );
      else
      {
        this->cycle();
        this->bcx.fill_halos(
          this->psi[this->n], ext(this->j, this->hlo)
        );
        this->bcy.fill_halos(
          this->psi[this->n], ext(this->i, this->hlo)
        );

        // choosing input/output for antidiff C
        const arrvec_t 
          &C_unco = (step == 1) 
            ? this->C 
            : (step % 2) 
              ? tmp[1]  // odd steps
              : tmp[0], // even steps
          &C_corr = (step  % 2) 
            ? tmp[0]    // odd steps
            : tmp[1];   // even steps

        // calculating the antidiffusive C 
        C_corr[0](im+h, this->j) = mpdata_C_adf<0>(
          this->psi[this->n], im, this->j, C_unco
        );
        this->bcy.fill_halos(C_corr[0], ext(this->i,h));

        C_corr[1](this->i, jm+h) = mpdata_C_adf<1>(
          this->psi[this->n], jm, this->i, C_unco
        );
        this->bcx.fill_halos(C_corr[1], ext(this->j,h));

        // donor-cell step 
        donorcell_op(
          this->psi, this->n, C_corr, this->i, this->j
        );
      }
    }
  }
};
//listing20
