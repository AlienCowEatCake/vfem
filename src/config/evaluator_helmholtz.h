#ifndef EVALUATOR_HELMHOLTZ_H
#define EVALUATOR_HELMHOLTZ_H

#include <complex>
#include <string>
#include "evaluator/evaluator.h"

// Вычислитель, заточенный под параметры уравнения Гельмгольца
class evaluator_helmholtz : public evaluator<std::complex<double> >
{
public:

    evaluator_helmholtz() : evaluator<std::complex<double> >()
    {
        update_cache();
    }

    evaluator_helmholtz(const std::string & str) : evaluator<std::complex<double> >(str)
    {
        update_cache();
    }

    evaluator_helmholtz(const evaluator_helmholtz & other) : evaluator<std::complex<double> >(other)
    {
        update_cache();
    }

    const evaluator_helmholtz & operator = (const evaluator_helmholtz & other)
    {
        evaluator<std::complex<double> >::copy_from_other(other);
        update_cache();
        return * this;
    }

    inline void set_x(const std::complex<double> & x)
    {
        * m_x = x;
    }

    inline void set_y(const std::complex<double> & y)
    {
        * m_y = y;
    }

    inline void set_z(const std::complex<double> & z)
    {
        * m_z = z;
    }

    inline void set_J0(const std::complex<double> & J0)
    {
        * m_J0 = J0;
    }

    inline void set_omega(const std::complex<double> & omega)
    {
        * m_omega = omega;
    }

    inline void set_epsilon(const std::complex<double> & epsilon)
    {
        * m_epsilon = epsilon;
    }

    inline void set_sigma(const std::complex<double> & sigma)
    {
        * m_sigma = sigma;
    }

    inline void set_mu(const std::complex<double> & mu)
    {
        * m_mu = mu;
    }

    inline void set_k2(const std::complex<double> & k2)
    {
        * m_k2 = k2;
    }

private:

    std::complex<double> * m_x, * m_y, * m_z;
    std::complex<double> * m_J0, * m_omega, * m_epsilon, * m_sigma, * m_mu, * m_k2;

    void update_cache()
    {
        m_x = evaluator<std::complex<double> >::m_variables["x"].pointer();
        m_y = evaluator<std::complex<double> >::m_variables["y"].pointer();
        m_z = evaluator<std::complex<double> >::m_variables["z"].pointer();
        m_J0 = evaluator<std::complex<double> >::m_variables["J0"].pointer();
        m_omega = evaluator<std::complex<double> >::m_variables["omega"].pointer();
        m_epsilon = evaluator<std::complex<double> >::m_variables["epsilon"].pointer();
        m_sigma = evaluator<std::complex<double> >::m_variables["sigma"].pointer();
        m_mu = evaluator<std::complex<double> >::m_variables["mu"].pointer();
        m_k2 = evaluator<std::complex<double> >::m_variables["k2"].pointer();
    }
};

#endif // EVALUATOR_HELMHOLTZ_H

