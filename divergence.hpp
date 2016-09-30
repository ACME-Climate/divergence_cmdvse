
#ifndef _DIVERGENCE_HPP_
#define _DIVERGENCE_HPP_

#include <array>

#if defined(__INTEL_COMPILER)
#define NOVECDEP _Pragma("ivdep")
#define ALWAYSVECTORIZE _Pragma("vector always")
#define ALIGN(vardec) __declspec(align) vardec
#define ALIGNTO(vardec, boundary) \
  __declspec(align(boundary)) vardec
#define RESTRICT restrict
#elif defined(__GNUG__)
#if(__GNUG__ == 4 && __GNUC_MINOR__ >= 9) || __GNUG__ > 4
#define NOVECDEP _Pragma("GCC ivdep")
#define ALWAYSVECTORIZE _Pragma("GCC vector always")
#define RESTRICT restrict
#else
#pragma message( \
    "G++ <4.9 Does not support vectorization pragmas")
#define NOVECDEP
#define ALWAYSVECTORIZE
#define RESTRICT
#endif

#define ALIGN(vardec) __attribute__((aligned)) vardec
#define ALIGNTO(vardec, boundary) \
  __attribute__((aligned(boundary))) vardec
#endif

constexpr const int dim = 2;

template <int np, typename real>
using real_matrix = std::array<std::array<std::array<std::array<real, np>, np>, dim>, dim>;

template <int np, typename real>
using real_vector = std::array<std::array<std::array<real, np>, np>, dim>;

template <int np, typename real>
using real_scalar = std::array<std::array<real, np>, np>;

template <int np, typename real>
struct element {
  real_scalar<np, real> metdet;
  real_matrix<np, real> Dinv;
  real_scalar<np, real> rmetdet;
};

template <int np, typename real>
struct derivative {
  real_scalar<np, real> Dvv;
};

extern "C" {
using real = double;
void divergence_sphere_fortran(
    const real_vector<4, real> &RESTRICT,
    const derivative<4, real> &RESTRICT,
    const element<4, real> &RESTRICT,
    real_scalar<4, real> &RESTRICT);
}

template <int np, typename real>
__attribute__((noinline)) void divergence_sphere(
    const real_vector<np, real> &RESTRICT v,
    const derivative<np, real> &RESTRICT deriv,
    const element<np, real> &RESTRICT elem,
    real_scalar<np, real> &RESTRICT div) {
  /* Computes the spherical divergence of v based on the
   * provided metric terms in elem and deriv
   * Returns the divergence in div
   */
  using rs = real_scalar<np, real>;
  using rv = real_vector<np, real>;

  using rfv = real[np * np * dim];
  using rfs = real[np * np];
  using rfm = real[np * np * 2 * 2];

  const RESTRICT rfv &vflat = *reinterpret_cast<const rfv*>(&v);
  const RESTRICT rfs &Dvv = *reinterpret_cast<const rfs*>(&deriv.Dvv);
  const RESTRICT rfs &metdet = *reinterpret_cast<const rfs*>(&elem.metdet);
  const RESTRICT rfm &Dinv = *reinterpret_cast<const rfm*>(&elem.Dinv);
  const RESTRICT rfs &rmetdet = *reinterpret_cast<const rfs*>(&elem.rmetdet);
  RESTRICT rfs &divflat = *reinterpret_cast<rfs *>(&div);
  #if 1
  std::cout.precision(17);
  std::cout << "\n\nv:\n";
  for(int i = 0; i < dim; i++) {
    for(int j = 0; j < np; j++) {
      for(int k = 0; k < np; k++) {
	std::cout << v[i][j][k] << "     ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }

  std::cout << "\n\nDvv:\n";
  for(int j = 0; j < np; j++) {
    for(int k = 0; k < np; k++) {
      std::cout << deriv.Dvv[j][k] << "     ";
    }
    std::cout << "\n";
  }

  std::cout << "\n\nmetdet:\n";
  for(int j = 0; j < np; j++) {
    for(int k = 0; k < np; k++) {
      std::cout << elem.metdet[j][k] << "     ";
    }
    std::cout << "\n";
  }

  std::cout << "\n\nDinvFlat:\n";
  for(int i = 0; i < np * np * dim * dim; i++) {
    std::cout << Dinv[i] << "   ";
  }
  std::cout << "\n\nDinv:\n";
  for(int j = 0; j < dim; j++) {
    for(int k = 0; k < dim; k++) {
      for(int l = 0; l < np; l++) {
	for(int m = 0; m < np; m++) {
	  std::cout << elem.Dinv[j][k][l][m] << "     ";
	}
	std::cout << "\n";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  
  std::cout << "rmetdet:\n";
  for(int j = 0; j < np; j++) {
    for(int k = 0; k < np; k++) {
      std::cout << elem.rmetdet[j][k] << "     ";
    }
    std::cout << "\n";
  }
  std::cout << "\n\n";
  #endif

  /* Convert to contra variant form and multiply by g */
  ALIGNTO(RESTRICT rv gv, 16);
  RESTRICT rfv &gvflat = *reinterpret_cast<rfv *>(&gv);

  #if 1
  for(int i = 0; i < np; i++) {
    for(int j = 0; j < np; j++) {
      for(int k = 0; k < dim; k++) {
        gv[j][i][k] = elem.metdet[j][i] *
                      (elem.Dinv[j][i][k][0] * v[j][i][0] +
                       elem.Dinv[j][i][k][1] * v[j][i][1]);
      }
    }
  }
  #else
  for(int j = 0; j < np * np; j++) {
    gvflat[j * dim] = metdet[j] *
      (Dinv[j * dim * dim] * vflat[j * dim] +
       Dinv[j * dim * dim + 1] * v[j / np][j % np][1]);


    gvflat[j * dim + 1] = metdet[j] *
      (Dinv[j * dim * dim + dim] * v[j / np][j % np][0] +
       Dinv[j * dim * dim + dim + 1] * v[j / np][j % np][1]);
  }
  #endif
  /* Compute d/dx and d/dy */
  ALIGNTO(RESTRICT rs vvtemp, 16);
  #if 1
  for(int l = 0; l < np; l++) {
    for(int j = 0; j < np; j++) {
      ALIGNTO(real dudx00, 16) = 0.0;
      ALIGNTO(real dvdy00, 16) = 0.0;
      for(int i = 0; i < np; i++) {
        dudx00 += deriv.Dvv[l][i] * gv[j][i][0];
        dvdy00 += deriv.Dvv[l][i] * gv[i][j][1];
      }
      div[j][l] = dudx00;
      vvtemp[l][j] = dvdy00;
    }
  }
  #endif
  constexpr const real rrearth = 1.5683814303638645E-7;

  #if 1
  for(int i = 0; i < np; i++) {
    for(int j = 0; j < np; j++) {
      div[i][j] = (div[i][j] + vvtemp[i][j]) *
                  (elem.rmetdet[i][j] * rrearth);
    }
  }
  #endif
}

#endif
