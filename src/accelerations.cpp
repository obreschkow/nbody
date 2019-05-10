#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List accelerations(NumericVector m, NumericMatrix x, NumericMatrix a, double G, double rsmoothsqr) {

  // initialize variables
  int n = m.size();
  NumericVector dminsqr(n);
  for (int i = 0; i<n; i++) dminsqr(i) = 1e99;

  // compute accelerations
  for (int i = 0; i<(n-1); i++) {
    for (int j = i+1; j<n; j++) {
      if ((m(i)>0)||(m(j)>0)) {
        double dx = x(i,0)-x(j,0); // faster than using a vector dx()
        double dy = x(i,1)-x(j,1);
        double dz = x(i,2)-x(j,2);
        double dsqr = fmax(rsmoothsqr,dx*dx+dy*dy+dz*dz);
        double q = G/pow(dsqr,1.5);
        a(i,0) -= m[j]*q*dx; // faster than using vector operation on a(_,i)
        a(i,1) -= m[j]*q*dy;
        a(i,2) -= m[j]*q*dz;
        a(j,0) += m[i]*q*dx;
        a(j,1) += m[i]*q*dy;
        a(j,2) += m[i]*q*dz;
        dminsqr(i) = fmin(dminsqr(i),dsqr);
        dminsqr(j) = fmin(dminsqr(j),dsqr);
      }
    }
  }

  // determine adaptive time step
  double zmin = 1e99;
  for (int i = 0; i<n; i++) {
    double z = dminsqr(i)/(a(i,0)*a(i,0)+a(i,1)*a(i,1)+a(i,2)*a(i,2));
    if (z<zmin) zmin = z;
  }
  double dtvar = pow(zmin,0.25);

  // output
  List ret;
  ret["a"] = a;
  ret["dtvar"] = dtvar;
  return ret;
}
