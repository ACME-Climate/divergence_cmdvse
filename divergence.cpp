
#include <sys/time.h>
#include <fstream>
#include <iostream>

#include "divergence.hpp"
#include "timer/timer.hpp"

template <int np, typename real>
void readVelocity(real_vector<np, real> &v,
                  std::istream *input) {
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < np; j++) {
      for(int k = 0; k < np; k++) {
        (*input) >> v[k][j][i];
      }
    }
  }
}

template <int np, typename real>
void readElement(element<np, real> &elem,
                 std::istream *input) {
  for(int i = 0; i < np; i++) {
    for(int j = 0; j < np; j++) {
      (*input) >> elem.metdet[i][j];
      elem.rmetdet[i][j] = 1 / elem.metdet[i][j];
    }
  }
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      for(int k = 0; k < np; k++) {
        for(int l = 0; l < np; l++) {
          (*input) >> elem.Dinv[l][k][i][j];
        }
      }
    }
  }
}

template <int np, typename real>
void readDerivative(derivative<np, real> &deriv,
                    std::istream *input) {
  for(int i = 0; i < np; i++) {
    for(int j = 0; j < np; j++) {
      (*input) >> deriv.Dvv[j][i];
    }
  }
}

template <int np, typename real>
void readDivergence(real_scalar<np, real> &divergence,
                    std::istream *input) {
  for(int i = 0; i < np; i++) {
    for(int j = 0; j < np; j++) {
      (*input) >> divergence[i][j];
    }
  }
}

constexpr const int DIMS = 2;

template <int np, typename real>
void compareDivergences(
    const real_vector<np, real> &v,
    const element<np, real> &elem,
    const derivative<np, real> &deriv,
    const real_scalar<np, real> &divergence_e,
    const int numtests) {
  Timer::Timer time_c;
  /* Initial run to prevent cache timing from affecting us
   */
  real_scalar<np, real> divergence_c;
  for(int i = 0; i < numtests; i++) {
    divergence_sphere<np, real>(v, deriv, elem,
                                divergence_c);
  }

  time_c.startTimer();
  for(int i = 0; i < numtests; i++) {
    divergence_sphere<np, real>(v, deriv, elem,
                                divergence_c);
  }
  time_c.stopTimer();

  Timer::Timer time_f;
  real_scalar<np, real> divergence_f;
  for(int i = 0; i < numtests; i++) {
    divergence_sphere_fortran(v, deriv, elem, divergence_f);
  }
  time_f.startTimer();
  for(int i = 0; i < numtests; i++) {
    divergence_sphere_fortran(v, deriv, elem, divergence_f);
  }
  time_f.stopTimer();
  std::cout << "Divergence Errors\n";
  for(int i = 0; i < np; i++) {
    for(int j = 0; j < np; j++) {
      std::cout << divergence_c[i][j] - divergence_e[i][j]
                << "    "
                << divergence_f[i][j] - divergence_e[i][j]
                << "\n";
    }
    std::cout << "\n";
  }

  std::cout << "C++ Time:\n"
            << time_c << "\n\nFortran Time:\n"
            << time_f << "\n";
}

int main(int argc, char **argv) {
  using real = double;
  constexpr const int NP = 4;
  real_vector<NP, real> v;
  element<NP, real> elem;
  derivative<NP, real> deriv;
  real_scalar<NP, real> divergence_e;
  {
    std::istream *input;
    if(argc > 1) {
      input = new std::ifstream(argv[1]);
    } else {
      input = &std::cin;
    }
    readVelocity<NP, real>(v, input);
    readElement(elem, input);
    readDerivative(deriv, input);
    readDivergence<NP, real>(divergence_e, input);
    if(argc > 1) {
      delete input;
    }
  }

  constexpr const int defNumTests = 1e5;
  const int numtests =
      (argc > 2) ? std::stoi(argv[2]) : defNumTests;
  compareDivergences<NP, real>(v, elem, deriv, divergence_e,
                               numtests);
  return 0;
}
