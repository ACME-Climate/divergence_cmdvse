
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

constexpr const int dim = 2;

template <int np, typename real>
using real_vector = real[np][np][dim];

template <int np, typename real>
using real_scalar = real[np][np];

template <int np, typename real>
struct element {
  real_scalar<np, real> metdet;
  real Dinv[np][np][2][2];
  real_scalar<np, real> rmetdet;
};

template <int np, typename real>
struct derivative {
  real_scalar<np, real> Dvv;
};

extern "C" {
using real = double;
constexpr const int np = 4;

__declspec(concurrency_safe(profitable),const)
void divergence_sphere_fortran(
    const real v[restrict np][np][dim],
    const real Dvv[restrict np][np],
    const real metdet[restrict np][np],
    const real Dinv[restrict np][np][dim][dim],
    const real rmetdet[restrict np][np],
    real div[restrict np][np]);
}

template <int np, typename real>
__attribute__((noinline,concurrency_safe(profitable))) void divergence_sphere(
    const real v[restrict np][np][dim],
    const real Dvv[restrict np][np],
    const real metdet[restrict np][np],
    const real Dinv[restrict np][np][dim][dim],
    const real rmetdet[restrict np][np],
    real div[restrict np][np]) {
  /* Computes the spherical divergence of v based on the
   * provided metric terms in metdet, Dinv, and Dvv
   * Returns the divergence in div
   */
  using rs = real_scalar<np, real>;
  using rv = real_vector<np, real>;
  constexpr const int alignment = 64;
  /* Convert to contra variant form and multiply by g */
  ALIGNTO(restrict rv gv, alignment);
  #if 1
  for(int j = 0; j < np; j++) {
    for(int i = 0; i < np; i++) {
      for(int k = 0; k < dim; k++) {
        gv[j][i][k] = metdet[j][i] *
                      (Dinv[j][i][k][0] * v[j][i][0] +
                       Dinv[j][i][k][1] * v[j][i][1]);
      }
    }
  }
  #endif
  /* Compute d/dx and d/dy */
  ALIGNTO(restrict rs vvtemp, alignment);
  #if 1
  for(int l = 0; l < np; l++) {
    for(int j = 0; j < np; j++) {
      ALIGNTO(real dudx00, alignment) = 0.0;
      ALIGNTO(real dvdy00, alignment) = 0.0;
      for(int i = 0; i < np; i++) {
        dudx00 = dudx00 + Dvv[l][i] * gv[j][i][0];
        dvdy00 += Dvv[l][i] * gv[i][j][1];
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
                  (rmetdet[i][j] * rrearth);
    }
  }
  #endif
}

#endif
