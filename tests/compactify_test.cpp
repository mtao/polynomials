#include "univariatepolynomial.h"
#include "univariatelegendre.h"
#include <cassert>


int main(int argc, char * argv[]) {
    UnivariatePolynomial<double> a(5);
    a(0) = 5;
    a(3) = 1;
    a.compactify();
     assert (a(0) == 5 );
     assert (a(1) == 0 );
     assert (a(2) == 0 );
     assert (a(3) == 1);
     assert (a.size() == 4);

     return 0;

}
