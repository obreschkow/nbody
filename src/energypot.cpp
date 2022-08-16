#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List energypot(NumericVector m, NumericMatrix x, NumericMatrix v, double rsmooth) {

  int n = m.size();
  NumericVector Epot(n);
  double dEpot;

  for (int i = 0; i<(n-1); i++) {
    for (int j = i+1; j<n; j++) {
      if ((m(i)>0)&&(m(j)>0)) {
        double dx = x(i,0)-x(j,0); // faster than using a vector dx()
        double dy = x(i,1)-x(j,1);
        double dz = x(i,2)-x(j,2);
        double r = sqrt(dx*dx+dy*dy+dz*dz);
        if (r>rsmooth) {
          dEpot = -m(i)*m(j)/r;
        } else {
          dEpot = -m(i)*m(j)*(3-pow(r/rsmooth,2))/(2*rsmooth);
        }
        Epot(i) = Epot(i)+dEpot;
        Epot(j) = Epot(j)+dEpot;
      }
    }
  }

  List ret;
  ret["Epot"] = Epot;
  return ret;
}
