
#ifndef _DIVERGENCE_HPP_
#define _DIVERGENCE_HPP_

#if defined(__INTEL_COMPILER)
#define NOVECDEP _Pragma("ivdep")
#define ALWAYSVECTORIZE _Pragma("vector always")
#define ALIGN(vardec) __declspec(align) vardec
#define ALIGNTO(vardec, boundary) \
  __declspec(align(boundary)) vardec
#elif defined(__GNUG__)
#if(__GNUG__ == 4 && __GNUC_MINOR__ >= 9) || __GNUG__ > 4
#define NOVECDEP _Pragma("GCC ivdep")
#define ALWAYSVECTORIZE _Pragma("GCC vector always")
#else
#pragma message( \
    "G++ <4.9 Does not support vectorization pragmas")
#define NOVECDEP
#define ALWAYSVECTORIZE
#endif

#define ALIGN(vardec) __attribute__((aligned)) vardec
#define ALIGNTO(vardec, boundary) \
  __attribute__((aligned(boundary))) vardec
#endif

template <int np, typename real>
struct element {
  real metdet[np][np];
  real Dinv[np][np][2][2];
  real rmetdet[np][np];
};

template <int np, typename real>
struct derivative {
  real Dvv[np][np];
};

constexpr const int dim = 2;

extern "C" {
using real = double;
void divergence_sphere_fortran(const real (*)[4][2],
                               const derivative<4, real> &,
                               const element<4, real> &,
                               real[4][4]);
}

template <int np, typename real>
__attribute__((noinline)) void divergence_sphere(
    const real v[np][np][dim],
    const derivative<np, real> &deriv,
    const element<np, real> &elem, real div[np][np]) {
  /* Convert to contra variant form and multiply by g */
  ALIGNTO(real gv[np][np][dim], 16);
  for(int j = 0; j < np; j++) {
    for(int i = 0; i < np; i++) {
      for(int k = 0; k < dim; k++) {
        gv[j][i][k] = elem.metdet[j][i] *
                      (elem.Dinv[j][i][k][0] * v[j][i][0] +
                       elem.Dinv[j][i][k][1] * v[j][i][1]);
      }
    }
  }
  /* Compute d/dx and d/dy */
  ALIGNTO(real vvtemp[np][np], 16);
  NOVECDEP
  for(int j = 0; j < np; j++) {
    NOVECDEP
    for(int l = 0; l < np; l++) {
      real dudx00 = 0.0;
      real dvdy00 = 0.0;
      NOVECDEP
      for(int i = 0; i < np; i++) {
        dudx00 += deriv.Dvv[l][i] * gv[j][i][0];
        dvdy00 += deriv.Dvv[l][i] * gv[i][j][1];
      }
      div[j][l] = dudx00;
      vvtemp[l][j] = dvdy00;
    }
  }
  constexpr const real rrearth = 1.5683814303638645E-7;

  for(int i = 0; i < np; i++) {
    for(int j = 0; j < np; j++) {
      div[i][j] = (div[i][j] + vvtemp[i][j]) *
                  (elem.rmetdet[i][j] * rrearth);
    }
  }
}

#endif
