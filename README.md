# Cubic Root Solver

This is part of the CPEN 333 lecture on using Git.

The goal of this challenge is to create a solver for finding all
the unique roots of a cubic polynomial.

 You are provided with a header-file that defines a few useful methods:
 ```cpp
 // solves for roots of a*x + b = 0,
 // returns number of roots found
 int solve_linear(double b, double a, double out[]);

 // solves for unique roots of a*x^2 + b*x + c = 0,
 // returns number of unique roots found
 int solve_quadratic(double c, double b, double a, double out[]);

 // solves for a root in the given range
 // requires f(left)*f(right) < 0, where
 //   f = sum(coeffs[i]*x^i, i=0..ncoeffs-1)
 // i.e. the root must be bracketed
 int solve_in_range(double coeffs[], int ncoeffs, double left,
                    double right, double &out)
 ```

There are several approaches for solving for roots of a cubic polynomial.
We will describe one of them here.

If we are able to partition the space into three regions,
  - `[-infinity, x1]`,
  - `[x1, x2]`,
  - `[x2, +infinity]`,

with `x1 <= x2`, and if we can guarantee that at most one root exists in
each partition, then all we would need to do is *check out each of the
branches for a possible solution*.

How do we find such an `x1` and `x2`?  There can be at most one root
in any region of a function where it is monotically increasing or
monotically decreasing.  These regions are always separated by what
are known as *stationary points*, where the derivative of the function
(in this case, the polynomial) is zero.  The derivative of a cubic
polynomial is a quadratic, so we should easily
be able to find the locations of the stationary points.

If there is a root in the monotonically increasing or monotically decreasing
region `[a, b]`, then we must have that `f(a)` and `f(b)` have opposite
signs.  Otherwise, since `f` is assumed monotone in the region, `f(x)` is
restricted to the region `[f(a), f(b)]` for `x` in `[a, b]`, and hence
would be bounded away from zero if they have the same sign.
If `f(a)` and `f(b)` have opposite signs, then
we can we can use the [bisection method](https://en.wikipedia.org/wiki/Bisection_method)
to find the root.

The solution approach is therefore as follows:

- Find the stationary points `x1` and `x2` of the cubic polynomial
- *Check out the branch* `[-infinity, x1]` for a possible solution using the
  bisection method (you may need to first find an appropriate lower-bound)
- *Check out the branch* `[x1, x2]` for a possible solution
- *Check out the branch* `[x2, infinity]` for a possible solution

The main program is provided to you.  Your job is to fill in the function
`int solve_cubic(...)`.