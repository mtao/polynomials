#include <univariatepolynomial.h>
#include <iostream>
#include <cassert>


int main(int argc, char * argv[]) {
    UnivariatePolynomial<double> a({1,2,3,4,5,6});
    assert(a.i() == UnivariatePolynomial<double>({0,1,1,1,1,1,1}));
    assert(a.d() == UnivariatePolynomial<double>({1*2,2*3,3*4,4*5,5*6}));

}
