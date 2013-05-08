#ifndef UNIVAR_POLY_H
#define UNIVAR_POLY_H
#include <Eigen/Dense>
#include <iomanip>

template <typename Scalar>
class UnivariatePolynomial {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        typedef typename Eigen::Matrix<Scalar,Eigen::Dynamic,1> Vector;
        UnivariatePolynomial(int size=0): m_coefficients(Vector::Zero(size+1)) {}
        UnivariatePolynomial(const Eigen::Ref<const Vector> & v): m_coefficients(v) {compactify();}
        UnivariatePolynomial(std::initializer_list<Scalar> l): m_coefficients(l.size()) {
            int i=0;
            for(const Scalar & v: l) {
                m_coefficients(i++) = v;
            }
            compactify();
        }


        void compactify() {
            int i=m_coefficients.rows()-1;
            if(i==0) return;
            while(i > 1 && m_coefficients(i) == Scalar(0))
            {
                --i;
            }
            m_coefficients = m_coefficients.head(i+1).eval();
        }

        UnivariatePolynomial<Scalar> d() {
            UnivariatePolynomial<Scalar> ret(Vector::Zero(size()-1));
            for(int i=1; i < this->size(); ++i) {
                ret(i-1) = i * this->coeffs()(i);
            }

        }
        Scalar eval(Scalar x) {
            Scalar ret=0;
            for(int i=coeffs().rows()-1; i >= 0; --i) {
                ret = x*(coeffs()(i) + ret);
            }
            return ret;
        }

        inline int size() const {
            return m_coefficients.rows();
        }
        inline int degree() const {
            return size()+1;
        }
        const Vector & coeffs() const {
            return m_coefficients;
        }
        Scalar operator()(int ind) const {return m_coefficients(ind);}
        Scalar & operator()(int ind) {return m_coefficients(ind);}
    private:
        Vector m_coefficients;
};

template <typename Scalar1, typename Scalar2>
UnivariatePolynomial<Scalar1> operator*(Scalar2 m, const UnivariatePolynomial<Scalar1> &a) {
    return UnivariatePolynomial<Scalar1>(m*a.coeffs());
}

template <typename Scalar>
UnivariatePolynomial<Scalar> operator*(const UnivariatePolynomial<Scalar> & lhs,const UnivariatePolynomial<Scalar> & rhs) {
    typedef typename UnivariatePolynomial<Scalar>::Vector Vector;
    int maxsize = lhs.size()+rhs.size();
    Vector v(Vector::Zero(maxsize));
    for(int i=0; i < lhs.size(); ++i) {
        for(int j=0; j < rhs.size(); ++j) {
            v(i+j) += lhs(i) * rhs(j);
        }
    }
    return UnivariatePolynomial<Scalar>(v);
}
template <typename Scalar>
UnivariatePolynomial<Scalar> operator+(const UnivariatePolynomial<Scalar> & lhs,const UnivariatePolynomial<Scalar> & rhs) {
    typedef typename UnivariatePolynomial<Scalar>::Vector Vector;
    int maxsize = std::max(lhs.size(),rhs.size());
    Vector v(Vector::Zero(maxsize));
    v.head(lhs.size()) = lhs.coeffs();
    v.head(rhs.size()).noalias() += rhs.coeffs();

    return UnivariatePolynomial<Scalar>(v);
}
template <typename Scalar>
UnivariatePolynomial<Scalar> operator-(const UnivariatePolynomial<Scalar> & lhs,const UnivariatePolynomial<Scalar> & rhs) {
    typedef typename UnivariatePolynomial<Scalar>::Vector Vector;
    int maxsize = std::max(lhs.size(),rhs.size());
    Vector v(Vector::Zero(maxsize));
    v.head(lhs.size()) = lhs.coeffs();
    v.head(rhs.size()).noalias() -= rhs.coeffs();

    return UnivariatePolynomial<Scalar>(v);
}
template <typename Scalar>
UnivariatePolynomial<Scalar> operator-(const UnivariatePolynomial<Scalar> & other) {
    return UnivariatePolynomial<Scalar>(-other->coeffs());
}

#endif
