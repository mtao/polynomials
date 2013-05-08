#include <univariatelegendre.h>
#include <iostream>
#include <cassert>


int main(int argc, char * argv[]) {
    UnivariatePolynomial<double> a({1,2,3,4,5,6});
    LegendrePolynomials<double> leg(10);
    assert((leg.getPoly(leg.project(a))-a).coeffs().norm() < 1e-10);
    return 0;


}
