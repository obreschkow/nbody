#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List accelerations(NumericVector m, NumericMatrix x, NumericMatrix v, NumericMatrix a, double G, double rsmoothsqr) {

  // initialize variables
  int n = m.size();
  NumericVector mindxsqr(n);
  for (int i = 0; i<n; i++) mindxsqr(i) = 1e99;
  //double w = 2e99;

  // iterate over point-pairs
  for (int i = 0; i<(n-1); i++) {
    for (int j = i+1; j<n; j++) {
      if ((m(i)>0)||(m(j)>0)) {

        // compute accelerations
        double dx0 = x(i,0)-x(j,0); // faster than using a vector dx0()
        double dx1 = x(i,1)-x(j,1);
        double dx2 = x(i,2)-x(j,2);
        double dxsqr = fmax(rsmoothsqr,dx0*dx0+dx1*dx1+dx2*dx2);
        double q = G/pow(dxsqr,1.5);
        a(i,0) -= m[j]*q*dx0; // faster than using vector operation on a(_,i)
        a(i,1) -= m[j]*q*dx1;
        a(i,2) -= m[j]*q*dx2;
        a(j,0) += m[i]*q*dx0;
        a(j,1) += m[i]*q*dx1;
        a(j,2) += m[i]*q*dx2;

        // partial time step computations
        mindxsqr(i) = fmin(mindxsqr(i),dxsqr);
        mindxsqr(j) = fmin(mindxsqr(j),dxsqr);
        //double dv0 = v(i,0)-v(j,0); // faster than using a vector dx0()
        //double dv1 = v(i,1)-v(j,1);
        //double dv2 = v(i,2)-v(j,2);
        //double dvsqr = dv0*dv0+dv1*dv1+dv2*dv2+1e-50;
        //w = fmin(w,dxsqr/dvsqr);

      }
    }
  }

  // finalize adaptive time step
  //double dtv = pow(w,0.5)/2; // recommended time step considering velocities
  double z = 1e99;
  for (int i = 0; i<n; i++) {
    double asqr = a(i,0)*a(i,0)+a(i,1)*a(i,1)+a(i,2)*a(i,2)+1e-50;
    z = fmin(z,mindxsqr(i)/asqr);
  }
  double dta = pow(z,0.25); // recommended time step considering accelerations
  double dtvar = dta;//fmin(1e99*dtv,dta);

  // output
  List ret;
  ret["a"] = a;
  ret["dtvar"] = dtvar;
  return ret;
}
