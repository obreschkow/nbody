#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List accelerations(NumericVector m, NumericMatrix x, NumericMatrix v, NumericMatrix a, double G, double rsmoothsqr, double L) {

  // initialize variables
  int n = m.size();
  NumericVector mindxsqr(n);
  for (int i = 0; i<n; i++) mindxsqr(i) = 1e99;
  double h = L/2.0;

  // iterate over point-pairs

  if (min(m)>0) {

    for (int i = 0; i<(n-1); i++) {
      for (int j = i+1; j<n; j++) {

        // compute accelerations
        double dx = x(i,0)-x(j,0); // faster than using a vector dx()
        double dy = x(i,1)-x(j,1);
        double dz = x(i,2)-x(j,2);
        if (L>0.0) {
          dx = std::fmod(dx+h+L,L)-h;
          dy = std::fmod(dy+h+L,L)-h;
          dz = std::fmod(dz+h+L,L)-h;
        }
        double dxsqr = fmax(rsmoothsqr,dx*dx+dy*dy+dz*dz);
        double q = G/pow(dxsqr,1.5);
        a(i,0) -= m[j]*q*dx; // faster than using vector operation on a(_,i)
        a(i,1) -= m[j]*q*dy;
        a(i,2) -= m[j]*q*dz;
        a(j,0) += m[i]*q*dx;
        a(j,1) += m[i]*q*dy;
        a(j,2) += m[i]*q*dz;

        // partial time step computations
        mindxsqr(i) = fmin(mindxsqr(i),dxsqr);
        mindxsqr(j) = fmin(mindxsqr(j),dxsqr);

      }
    }

  } else {

    // iterate over point-pairs
    for (int i = 0; i<n; i++) { // acceleration of particle i
      if (m(i)>=0) {
        for (int j = 0; j<n; j++) { // ... due to particle j
          if (j!=i) {
            double mj = fabs(m(j));
            if (mj>0) {

              // compute accelerations
              double dx = x(j,0)-x(i,0); // faster than using a vector dx()
              double dy = x(j,1)-x(i,1);
              double dz = x(j,2)-x(i,2);
              if (L>0.0) {
                dx = std::fmod(dx+h+L,L)-h;
                dy = std::fmod(dy+h+L,L)-h;
                dz = std::fmod(dz+h+L,L)-h;
              }
              double dxsqr = fmax(rsmoothsqr,dx*dx+dy*dy+dz*dz);
              double q = G/pow(dxsqr,1.5);
              a(i,0) += mj*q*dx; // faster than using vector operation on a(_,i)
              a(i,1) += mj*q*dy;
              a(i,2) += mj*q*dz;

              // partial time step computations
              mindxsqr(i) = fmin(mindxsqr(i),dxsqr);

            }
          }
        }
      }
    }
  }

  // finalize adaptive time step
  double z = 1e99;
  for (int i = 0; i<n; i++) {
    double asqr = a(i,0)*a(i,0)+a(i,1)*a(i,1)+a(i,2)*a(i,2)+1e-50;
    z = fmin(z,mindxsqr(i)/asqr);
  }
  double dtvar = pow(z,0.25); // recommended time step considering accelerations

  // output
  List ret;
  ret["a"] = a;
  ret["dtvar"] = dtvar;
  return ret;
}
