#include "univariatepolynomial.h"

#include <vector>
#include <iostream>
#include <cassert>

template <typename Scalar>
class LegendrePolynomials {
    public:
    typedef UnivariatePolynomial<Scalar> Polynomial;
    typedef typename Polynomial::Vector Vector;
    LegendrePolynomials(int depth=2): m_polynomials(depth) {
            assert(depth > 2);

        m_polynomials[0]=Polynomial({1.0});
        m_polynomials[1]=Polynomial({0,1});
        Polynomial x({0,1});
        for(int n=2; n < depth; ++n) {
            m_polynomials[n] = 1.0/Scalar(n) * ((2*n-1) * x * m_polynomials[n-1] - (n-1) * m_polynomials[n-2]);
            Eigen::Matrix<Scalar,Eigen::Dynamic,1> m = m_polynomials[n].coeffs();
        }
    }
    const Polynomial & get(int depth){return m_polynomials[depth];}

    Vector project(Polynomial f) {
        int maxsize = m_polynomials.size();
        Vector v(maxsize);
        for(int i=0; i < maxsize; ++i) {
            auto& b = m_polynomials[i];
            v(i) = (f*b).i(-1,1) / (b*b).i(-1,1);
        }
        return v;
    }

    Polynomial getPoly(const Vector & coeffs) {
        Polynomial p;
        for(int i=0; i < m_polynomials.size(); ++i) {
            p += coeffs(i) * m_polynomials[i];
            
        }
        return p;
    }



    private:
    std::vector<Polynomial> m_polynomials;
};
