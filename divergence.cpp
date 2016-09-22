
#include <iostream>
#include <fstream>
#include <sys/time.h>

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
extern void divergence_sphere_fortran(const double (*)[4][2],
                                      const derivative<4, double> &,
                                      const element<4, double> &, double[4][4]);
}

template <int np, typename real>
void divergence_sphere(const real v[np][np][2],
                       const derivative<np, real> &deriv,
                       const element<np, real> &elem, real div[np][np]) {
  real gv[np][np][2];
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < np; j++) {
      gv[i][j][0] = elem.metdet[i][j] * (elem.Dinv[i][j][0][0] * v[i][j][0] +
                                         elem.Dinv[i][j][0][1] * v[i][j][1]);
      gv[i][j][1] = elem.metdet[i][j] * (elem.Dinv[i][j][1][0] * v[i][j][0] +
                                         elem.Dinv[i][j][1][1] * v[i][j][1]);
    }
  }
  real vvtemp[np][np];
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < np; j++) {
      real dudx00 = 0.0;
      real dvdy00 = 0.0;
      for (int k = 0; k < np; k++) {
        dudx00 += deriv.Dvv[j][k] * gv[i][j][0];
        dvdy00 += deriv.Dvv[j][k] * gv[i][j][1];
      }
      div[i][j] = dudx00;
      vvtemp[i][j] = dvdy00;
    }
  }
  constexpr const real rrearth = 1.5683814303638645E-7;
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < np; j++) {
      div[i][j] = (div[i][j] + vvtemp[i][j]) * (elem.rmetdet[i][j] * rrearth);
    }
  }
}

template <int np, typename real>
void readVelocity(real v[np][np][2], std::istream *input) {
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < np; j++) {
      for (int k = 0; k < np; k++) {
        (*input) >> v[i][j][k];
      }
    }
  }
}

template <int np, typename real>
void readElement(element<np, real> &elem, std::istream *input) {
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < np; j++) {
      (*input) >> elem.metdet[i][j];
    }
  }
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < np; j++) {
      (*input) >> elem.rmetdet[j][i];
    }
  }
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < np; k++) {
        for (int l = 0; l < np; l++) {
          (*input) >> elem.Dinv[k][l][i][j];
        }
      }
    }
  }
}

template <int np, typename real>
void readDerivative(derivative<np, real> &deriv, std::istream *input) {
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < np; j++) {
      (*input) >> deriv.Dvv[j][i];
    }
  }
}

int main(int argc, char **argv) {
  std::istream *input;
  if (argc > 1) {
    input = new std::ifstream(argv[1]);
  } else {
    input = &std::cin;
  }
  using real = double;
  constexpr const int DIMS = 2;
  constexpr const int NP = 4;
  double v[NP][NP][DIMS];
  readVelocity(v, input);
  struct element<NP, real> elem;
  readElement(elem, input);
  struct derivative<NP, real> deriv;
  readDerivative(deriv, input);
  constexpr const char *dimname[] = {"U", "V"};
  std::cout.precision(20);
  for (int dim = 0; dim < DIMS; dim++) {
    std::cout << "\n" << dimname[dim] << "\n";
    for (int i = 0; i < NP; i++) {
      for (int j = 0; j < NP; j++) {
        std::cout << "   " << v[i][j][dim] << "\n";
      }
    }
  }
  std::cout << "\nelem.metdet\n";
  for (int i = 0; i < NP; i++) {
    for (int j = 0; j < NP; j++) {
      std::cout << elem.metdet[i][j] << "\n";
    }
  }
  std::cout << "\nelem.rmetdet\n";
  for (int i = 0; i < NP; i++) {
    for (int j = 0; j < NP; j++) {
      std::cout << elem.rmetdet[i][j] << "\n";
    }
  }
  for (int i = 0; i < DIMS; i++) {
    for (int j = 0; j < DIMS; j++) {
      std::cout << "\nelem.Dinv (" << i << ", " << j << ")\n";
      for (int k = 0; k < NP; k++) {
        for (int l = 0; l < NP; l++) {
          std::cout << elem.Dinv[k][l][i][j] << "\n";
        }
      }
    }
  }
  std::cout << "\nderiv.Dvv\n";
  for (int i = 0; i < NP; i++) {
    for (int j = 0; j < NP; j++) {
      std::cout << deriv.Dvv[i][j] << "\n";
    }
  }

  real divergence_c[NP][NP];
  divergence_sphere<NP, real>(v, deriv, elem, divergence_c);
  real divergence_f[NP][NP];
  divergence_sphere_fortran(v, deriv, elem, divergence_f);
  std::cout << "\nDivergence\n";
  for (int i = 0; i < NP; i++) {
    for (int j = 0; j < NP; j++) {
      std::cout << divergence_c[i][j] << "    " << divergence_f[i][j] << "\n";
    }
  }
  if (argc > 1) {
    delete input;
  }
  return 0;
}
