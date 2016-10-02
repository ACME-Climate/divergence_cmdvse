
#ifndef _DIVERGENCE_HPP_
#define _DIVERGENCE_HPP_

#if defined(__INTEL_COMPILER)
#pragma message("Building with ICC")
#define NOVECDEP _Pragma("ivdep")
#define ALWAYSVECTORIZE _Pragma("vector always")
#define ALIGN(vardec) __declspec(align) vardec
#define ALIGNTO(vardec, boundary) \
  __declspec(align(boundary)) vardec
#define RESTRICT

#elif defined(__clang__)
#pragma message("Building with clang++")
#define NOVECDEP _Pragma("GCC ivdep")
#define ALWAYSVECTORIZE _Pragma("GCC vector always")
#define ALIGN(vardec) __attribute__((aligned)) vardec
#define ALIGNTO(vardec, boundary) \
  __attribute__((aligned(boundary))) vardec
#define RESTRICT

#elif defined(__GNUG__)
#pragma message("Building with g++")
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
#define RESTRICT __restrict__

#else
#pragma message("Unknown compiler")
#define ALIGN(vardec) vardec
#define ALIGNTO(vardec, boundary) vardec
#define NOVECDEP
#define ALWAYSVECTORIZE
#define RESTRICT
#endif

#include <array>

constexpr const int dim = 2;

template <int np, typename real>
using real_matrix = std::array<
    std::array<std::array<std::array<real, dim>, dim>, np>,
    np>;

template <int np, typename real>
using real_matrix_flat =
    std::array<real, dim * dim * np * np>;

template <int np, typename real>
using real_vector =
    std::array<std::array<std::array<real, dim>, np>, np>;

template <int np, typename real>
using real_vector_flat = std::array<real, dim * np * np>;

template <int np, typename real>
using real_scalar = std::array<std::array<real, np>, np>;

template <int np, typename real>
using real_scalar_flat = std::array<real, np * np>;

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
    const real_vector<np, real> RESTRICT &v,
    const derivative<np, real> &RESTRICT deriv,
    const element<np, real> &RESTRICT elem,
    real_scalar<np, real> RESTRICT &div) {
#if 0
  std::cout.precision(17);
	std::cout.width(20);
  std::cout << "\nv:\n";
  for(int i = 0; i < np; i++) {
    for(int j = 0; j < np; j++) {
      for(int k = 0; k < dim; k++) {
				std::cout << v[i][j][k] << "  ";
      }
			std::cout << "\n";
    }
		std::cout << "\n";
  }

  std::cout << "\nDvv:\n";
  for(int j = 0; j < np; j++) {
    for(int k = 0; k < np; k++) {
      std::cout << deriv.Dvv[j][k] << "  ";
    }
		std::cout << "\n";
  }

  std::cout << "\nmetdet:\n";
  for(int j = 0; j < np; j++) {
    for(int k = 0; k < np; k++) {
      std::cout << elem.metdet[j][k] << "  ";
    }
		std::cout << "\n";
  }

  std::cout << "\nDinv:\n";
  for(int j = 0; j < np; j++) {
    for(int k = 0; k < np; k++) {
      for(int l = 0; l < dim; l++) {
				for(int m = 0; m < dim; m++) {
					std::cout << elem.Dinv[j][k][l][m] << "     ";
				}
				std::cout << "\n";
      }
			std::cout << "\n";
    }
		std::cout << "\n";
  }
#endif
  /* Computes the spherical divergence of v based on the
   * provided metric terms in elem and deriv
   * Returns the divergence in div
   */
  using rs = real_scalar<np, real>;
  using rsf = real_scalar_flat<np, real>;
  using rv = real_vector<np, real>;
  using rvf = real_vector_flat<np, real>;
  /* Convert to contra variant form and multiply by g */
  ALIGNTO(rvf RESTRICT gv, 16);
#if 1
  for(int j = 0; j < np; j++) {
    for(int i = 0; i < np; i++) {
      for(int k = 0; k < dim; k++) {
        gv[(j * np + i) * dim + k] =
            elem.metdet[j][i] *
            (elem.Dinv[j][i][k][0] * v[j][i][0] +
             elem.Dinv[j][i][k][1] * v[j][i][1]);
      }
    }
  }
#endif
  /* Compute d/dx and d/dy */
  ALIGNTO(rs RESTRICT vvtemp, 16);
#if 1
  for(int l = 0; l < np; l++) {
    for(int j = 0; j < np; j++) {
      ALIGNTO(real dudx00, 16) = 0.0;
      ALIGNTO(real dvdy00, 16) = 0.0;
      for(int i = 0; i < np; i++) {
        dudx00 +=
            deriv.Dvv[l][i] * gv[(j * np + i) * dim + 0];
        dvdy00 +=
            deriv.Dvv[l][i] * gv[(i * np + j) * dim + 1];
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
