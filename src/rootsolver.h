#ifndef ROOTSOLVER_H
#define ROOTSOLVER_H

#include <cmath>

namespace rootsolver {

// solves for roots of a*x + b = 0,
// returns number of roots found
int solve_linear(double b, double a, double out[]) {  
  if (a == 0) {
    return 0;
  }
  out[0] = -b/a;
  return 1;
}

// solves for unique roots of a*x^2 + b*x + c = 0,
// returns number of unique roots found
int solve_quadratic(double c, double b, double a, double out[]) {

  // check if actually linear
  if (a == 0) {
    return solve_linear(b, c, out);
  }

  double q = b*b-4*a*c;

  // imaginary
  if (q < 0) {
    return 0;
  }
  
  // single real root
  if (q == 0) {
    // -b/2a   
    out[0] = -b/(2*a);
    return 1;
  }

  // two real roots
  double rootq = std::sqrt(q);
  out[0] = (-b-rootq)/(2*a);
  out[1] = (-b+rootq)/(2*a);
  return 2;
}

// solves for a root in the given range
// requires f(left)*f(right) < 0, where f = sum(coeffs[i]*x^i, i=0..ncoeffs-1)
//   i.e. the root must be bracketed
int solve_in_range(double coeffs[], int ncoeffs, double left, double right, double &out) {

  // binary search in range

  // compute f(left)
  double fleft = 0;
  for (int i=0; i<ncoeffs; ++i) {
    fleft += coeffs[i]*std::pow(left, i);
  }

  // compute f(right)
  double fright = 0;
  for (int i=0; i<ncoeffs; ++i) {
    fright += coeffs[i]*std::pow(right, i);
  }

  // root not necessarily bracketed
  if (fleft*fright > 0) {
    return 0;
  } else if (fleft*fright == 0) {
    // return smaller of two
    if (std::fabs(fleft) < std::fabs(fright)) {
      return fleft;
    }
    return fright;
  }

  double mid = left;

  // accuracy of (right-left)/2^64
  for (int i=0; i<64; ++i) {
    mid = (left+right)/2;
    // compute f(mid)
    double fmid = 0;
    for (int i=0; i<ncoeffs; ++i) {
      fmid += coeffs[i]*std::pow(mid, i);
    }

    // check if we bracket on left or right
    if (fleft*fmid < 0) {
      right = mid;
      fright = fmid;
    } else if (fmid*fright < 0) {
      // we do this to enforce that fleft*fright < 0
      //   even with numerical oddities
      left = mid;
      fleft = fmid;
    } else {
      // funky area in-between, we've found root
      out = mid;
      return 1;
    }
  };

  out = mid;
  return 1;

}

} // rootsolver

#endif
