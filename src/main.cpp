#include "rootsolver.h"

#include <iostream>

// Solves for all unique roots of f(x) = a*x^3 + b*x^2 + c*x + d
//   populates out[] with the solution
//   returns # roots found
int solve_cubic(double d, double c, double b, double a, double out[]) {

   // YOUR CODE HERE //
   // See README.md for hints //
   
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
    << c << "x + " << d << " has " << nroots << " roots:"  ;
  for (int i=0; i<nroots; ++i) {
    std::cout << "  " << roots[i] << std::endl;
  }

  return 0;
}
