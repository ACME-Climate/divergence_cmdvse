
#include <iostream>
#include <fstream>
#include <sys/time.h>

#include "timer/timer.hpp"

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

extern "C" {
void divergence_sphere_fortran(
    const double (*)[4][2], const derivative<4, double> &,
    const element<4, double> &, double[4][4]);
}

template <int np, typename real>
void divergence_sphere(const real v[np][np][2],
                       const derivative<np, real> &deriv,
                       const element<np, real> &elem,
                       real div[np][np]) {
  /* Convert to contra variant form and multiply by g */
  real gv[np][np][2];
  for(int i = 0; i < np; i++) {
    for(int j = 0; j < np; j++) {
      gv[i][j][0] = elem.metdet[i][j] *
                    (elem.Dinv[i][j][0][0] * v[i][j][0] +
                     elem.Dinv[i][j][0][1] * v[i][j][1]);
      gv[i][j][1] = elem.metdet[i][j] *
                    (elem.Dinv[i][j][1][0] * v[i][j][0] +
                     elem.Dinv[i][j][1][1] * v[i][j][1]);
    }
  }

  /* Compute d/dx and d/dy */
  real vvtemp[np][np];
  for(int i = 0; i < np; i++) {
    for(int j = 0; j < np; j++) {
      real dudx00 = 0.0;
      real dvdy00 = 0.0;
      for(int k = 0; k < np; k++) {
        dudx00 += deriv.Dvv[j][k] * gv[i][k][0];
        dvdy00 += deriv.Dvv[j][k] * gv[k][i][1];
      }
      div[i][j] = dudx00;
      vvtemp[j][i] = dvdy00;
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

template <int np, typename real>
void readVelocity(real v[np][np][2], std::istream *input) {
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
    }
  }
  for(int i = 0; i < np; i++) {
    for(int j = 0; j < np; j++) {
      (*input) >> elem.rmetdet[j][i];
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

int main(int argc, char **argv) {
  using real = double;
  constexpr const int DIMS = 2;
  constexpr const int NP = 4;
  double v[NP][NP][DIMS];
  struct element<NP, real> elem;
  struct derivative<NP, real> deriv;
  {
    std::istream *input;
    if(argc > 1) {
      input = new std::ifstream(argv[1]);
    } else {
      input = &std::cin;
    }
    readVelocity(v, input);
    readElement(elem, input);
    readDerivative(deriv, input);
    if(argc > 1) {
      delete input;
    }
  }
  
  constexpr const int numtests = 1e5;
  Timer::Timer time_c;
  /* Initial run to prevent cache timing from affecting us */
  real divergence_c[NP][NP];
  for(int i = 0; i < numtests; i++) {
    divergence_sphere<NP, real>(v, deriv, elem, divergence_c);
  }
  
  time_c.startTimer();
  for(int i = 0; i < numtests; i++) {
    divergence_sphere<NP, real>(v, deriv, elem, divergence_c);
  }
  time_c.stopTimer();
  
  Timer::Timer time_f;
  real divergence_f[NP][NP];
  for(int i = 0; i < numtests; i++) {
    divergence_sphere_fortran(v, deriv, elem, divergence_f);
  }
  time_f.startTimer();
  for(int i = 0; i < numtests; i++) {
    divergence_sphere_fortran(v, deriv, elem, divergence_f);
  }
  time_f.stopTimer();
  constexpr const real divergence_e[NP][NP] = {
      {0.14383368343270220, 0.16634973900122296,
       0.21556384471655873, 0.29459485031864308},

      {0.10495094678642192, 0.10563956451788356,
       0.13337591314546268, 0.18918612996975015},

      {0.10396974364110805, 0.10785355163561680,
       0.13773336288220531, 0.19410009772740583},

      {0.14264597462181144, 0.17154723246534875,
       0.22337488503182140, 0.30152141985185271}};
  std::cout << "Divergence Errors\n";
  for(int i = 0; i < NP; i++) {
    for(int j = 0; j < NP; j++) {
      std::cout << divergence_c[i][j] - divergence_e[i][j]
                << "    "
                << divergence_f[i][j] - divergence_e[i][j]
                << "\n";
    }
    std::cout << "\n";
  }

  std::cout << "C++ Time:\n" << time_c
	    << "\n\nFortran Time:\n" << time_f << "\n";
  return 0;
}
