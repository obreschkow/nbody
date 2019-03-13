#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List accelerations(NumericVector m, NumericMatrix x, double G, double rsmoothsqr) {

  // initialize variables
  int n = m.size();
  NumericMatrix a(3,n);
  NumericVector dminsqr(n);
  for (int i = 0; i<n; i++) dminsqr(i) = 1e99;

  // compute accelerations
  for (int i = 0; i<(n-1); i++) {
    for (int j = i+1; j<n; j++) {
      if ((m(i)>0)||(m(j)>0)) {
        double dx = x(0,i)-x(0,j); // faster than using a vector dx()
        double dy = x(1,i)-x(1,j);
        double dz = x(2,i)-x(2,j);
        double dsqr = fmax(rsmoothsqr,dx*dx+dy*dy+dz*dz);
        double q = G/pow(dsqr,1.5);
        a(0,i) -= m[j]*q*dx; // faster than using vector operation on a(_,i)
        a(1,i) -= m[j]*q*dy;
        a(2,i) -= m[j]*q*dz;
        a(0,j) += m[i]*q*dx;
        a(1,j) += m[i]*q*dy;
        a(2,j) += m[i]*q*dz;
        dminsqr(i) = fmin(dminsqr(i),dsqr);
        dminsqr(j) = fmin(dminsqr(j),dsqr);
      }
    }
  }

  // determine adaptive time step
  double zmin = 1e99;
  for (int i = 0; i<n; i++) {
    double z = dminsqr(i)/(a(0,i)*a(0,i)+a(1,i)*a(1,i)+a(2,i)*a(2,i));
    if (z<zmin) zmin = z;
  }
  double dtvar = pow(zmin,0.25);

  // output
  List ret;
  ret["a"] = a;
  ret["dtvar"] = dtvar;
  return ret;
}
