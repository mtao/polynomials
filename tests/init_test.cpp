#include <iostream>
#include "univariatepolynomial.h"
#include "univariatelegendre.h"


int main(int argc, char * argv[]) {
    UnivariatePolynomial<double> a(3);
    a(0) = 1;
    UnivariatePolynomial<double> b(3);
    b(1) = 1;
    UnivariatePolynomial<double> c({1,2,3,4,5,6,7});
    c(2) = 2;
    auto&& d = a + b + c;
    std::cout << a.eval(.5) << " " << b.eval(.5) << " " << c.eval(.5) << std::endl;
    std::cout << d.eval(.5) << std::endl;

    b.compactify();
    std::cout << "B size " << b.size() << std::endl;
    std::cout << (b * b).coeffs().transpose() << std::endl;

    std::cout << "Legendre: " << std::endl;

    LegendrePolynomials<double> lp(5);
    for(int i=0; i < 5; ++i) {
        std::cout << i << ") " << lp.get(i).coeffs().transpose() << std::endl;
    }
    return 0;
}
