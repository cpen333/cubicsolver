#include "rootsolver.h"

#include <iostream>

// Solves for all unique roots of f(x) = a*x^3 + b*x^2 + c*x + d
//   populates out[] with the solution
//   returns # roots found
int solve_cubic(double d, double c, double b, double a, double out[]) {

  // check if quadratic
  if (a == 0) {
    return rootsolver::solve_quadratic(d, c, b, out);
  }

  // roots of f(x) and -f(x) are the same, so for convenience enforce a>0
  if (a < 0) {
    a = -a;
    b = -b;
    c = -c;
    d = -d;
  }

  // store coefficients in array for solve_in_range
  double coeffs[4];
  coeffs[0] = d;
  coeffs[1] = c;
  coeffs[2] = b;
  coeffs[3] = a;

  // find stationary points
  int nstats = rootsolver::solve_quadratic(c, 2*b, 3*a, out);

  if (nstats == 0) {
    // no stationary points, go left until f<0 and right until f>0

    // find lower-bound
    double dx = -1;       // start with a certain step size to search left
    double x0 = dx;
    double f0 = a*x0*x0*x0 + b*x0*x0 + c*x0 + d;
    while (f0 > 0) {
      // keep doubling step size and move further left
      dx *= 2;
      x0 = x0+dx;
      f0 = a*x0*x0*x0 + b*x0*x0 + c*x0 + d;
    }

    // find upper-bound
    dx = 1;       // start with a certain step size to search right
    double x3 = dx;
    double f3 = a*x3*x3*x3 + b*x3*x3 + c*x3 + d;
    while (f3 < 0) {
      // keep doubling step size and move further right
      dx *= 2;
      x3 = x3+dx;
      f3 = a*x3*x3*x3 + b*x3*x3 + c*x3 + d;
    }

    return rootsolver::solve_in_range(coeffs, 4, x0, x3, out[0]);

  } else if (nstats == 1) {
    // inflection point, see if it's greater or less than zero
    double x1 = out[0];
    double f1 = a*x1*x1*x1 + b*x1*x1 + c*x1 + d;
    
    if (f1 == 0) {
      return 1;  // root at inflection point
    } else if (f1 < 0) {
      // single solution must be in [x1, infinity]
      // find upper-bound
      double dx = 1;       // start with a certain step size to search right
      double x3 = x1+dx;
      double f3 = a*x3*x3*x3 + b*x3*x3 + c*x3 + d;
      while (f3 < 0) {
        // keep doubling step size and move further right
        dx *= 2;
        x3 = x3+dx;
        f3 = a*x3*x3*x3 + b*x3*x3 + c*x3 + d;
      }

      // should have root in [x1, x3] now
      return rootsolver::solve_in_range(coeffs, 4, x1, x3, out[0]);
    } else {
      // single solution must be in [-infinity, x1]
      // find lower-bound
      double dx = -1;       // start with a certain step size to search left
      double x0 = x1+dx;
      double f0 = a*x0*x0*x0 + b*x0*x0 + c*x0 + d;
      while (f0 > 0) {
        // keep doubling step size and move further left
        dx *= 2;
        x0 = x0+dx;
        f0 = a*x0*x0*x0 + b*x0*x0 + c*x0 + d;
      }

      // should have root in [x0, x1] now
      return rootsolver::solve_in_range(coeffs, 4, x0, x1, out[0]);
    }
  } // single stationary point

  // we have two stationary points
  double x1 = out[0];
  double f1 = a*x1*x1*x1 + b*x1*x1 + c*x1 + d;
  double x2 = out[1];
  double f2 = a*x2*x2*x2 + b*x2*x2 + c*x2 + d;

  // left-most root between (-infinity, x1]
  int nroots = 0;
  if (f1 > 0) {
    // find lower-bound
    double dx = -1;       // start with a certain step size to search left
    double x0 = x1+dx;
    double f0 = a*x0*x0*x0 + b*x0*x0 + c*x0 + d;
    while (f0 > 0) {
      // keep doubling step size and move further left
      dx *= 2;
      x0 = x0+dx;
      f0 = a*x0*x0*x0 + b*x0*x0 + c*x0 + d;
    }

    nroots += rootsolver::solve_in_range(coeffs, 4, x0, x1, out[nroots]);
  } else if (f1 == 0) {
    // root exactly on x1
    out[nroots] = x1;
    x1 += 1e-16*(x2-x1);  // perturb root ever-so-slightly
    ++nroots;
  }

  // next root between (x1, x2)
  if (f1 > 0 && f2 < 0) {
    nroots += rootsolver::solve_in_range(coeffs, 4, x1, x2, out[nroots]);
  } 

  // right-most root between [x2, infinity)
  if (f2 == 0) {
    // root exactly at x2
    out[nroots] = x2;
    x2 += 1e-16*(x2-x1);  // perturb root ever-so-slightly
    ++nroots;
  } else if (f2 < 0) {
    // find upper-bound
    double dx = 1;       // start with a certain step size to search right
    double x3 = x2+dx;
    double f3 = a*x3*x3*x3 + b*x3*x3 + c*x3 + d;
    while (f3 < 0) {
      // keep doubling step size and move further right
      dx *= 2;
      x3 = x3+dx;
      f3 = a*x3*x3*x3 + b*x3*x3 + c*x3 + d;
    }

    nroots += rootsolver::solve_in_range(coeffs, 4, x2, x3, out[nroots]);
  }

  return nroots;

}

int main() {

  std::cout << "Cubic root solver" << std::endl;
  double a, b, c, d;
  std::cout << "  Enter the coefficient of x^3: ";
  std::cin >> a;
  std::cout << "  Enter the coefficient of x^2: ";
  std::cin >> b;
  std::cout << "  Enter the coefficient of x^1: ";
  std::cin >> c;
  std::cout << "  Enter the coefficient of x^0: ";
  std::cin >> d;
  std::cout << std::endl;

  double roots[3];  // storage for roots

  int nroots = solve_cubic(d, c, b, a, roots);
  std::cout << std::endl << a << "x^3 + " << b << "x^2 + " 
    << c << "x + " << d << " has " << nroots << " roots:" << std::endl;
  for (int i=0; i<nroots; ++i) {
    std::cout << "  " << roots[i] << std::endl;
  }

  return 0;
}
