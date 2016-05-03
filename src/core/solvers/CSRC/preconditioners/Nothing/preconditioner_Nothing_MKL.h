#if !defined(SOLVERS_CSRC_PRECONDITIONERS_NOTHING_MKL_H_INCLUDED)
#define SOLVERS_CSRC_PRECONDITIONERS_NOTHING_MKL_H_INCLUDED

#include <complex>
#include "preconditioner_Nothing.h"
#include "../../../../wrappers/omp_wrapper.h"
#include "../../../../wrappers/mkl_wrapper.h"

namespace core { namespace solvers { namespace CSRC { namespace preconditioners {

/**
 * @brief Класс базовый (пустой) предобуславливатель A = S * Q, использующий библиотеку MKL
 */
template<typename val_type, typename ind_type = std::size_t>
class preconditioner_Nothing_MKL : public preconditioner_Nothing<val_type, ind_type>
{
public:
    /**
     * @brief Конструктор для несимметричного предобуславливателя
     * @param[in] gi Массив ig в разряженном строчно-столбцовом представлении
     * @param[in] gj Массив jg в разряженном строчно-столбцовом представлении
     * @param[in] di Массив диагональных эл-тов в разряженном строчно-столбцовом представлении
     * @param[in] gl Массив нижнего треугольника в разряженном строчно-столбцовом представлении
     * @param[in] gu Массив верхнего треугольника в разряженном строчно-столбцовом представлении
     * @param[in] n Размерность матрицы
     */
    preconditioner_Nothing_MKL(const ind_type * gi, const ind_type * gj, const val_type * di,
                               const val_type * gl, const val_type * gu, ind_type n)
        : preconditioner_Nothing<val_type, ind_type>(gi, gj, di, gl, gu, n)
    {
        init();
    }

    /**
     * @brief Конструктор для симметричного предобуславливателя
     * @param[in] gi Массив ig в разряженном строчно-столбцовом представлении
     * @param[in] gj Массив jg в разряженном строчно-столбцовом представлении
     * @param[in] di Массив диагональных эл-тов в разряженном строчно-столбцовом представлении
     * @param[in] gg Массив нижнего и верхнего треугольника в разряженном строчно-столбцовом представлении
     * @param[in] n Размерность матрицы
     */
    preconditioner_Nothing_MKL(const ind_type * gi, const ind_type * gj, const val_type * di,
                               const val_type * gg, ind_type n)
        : preconditioner_Nothing<val_type, ind_type>(gi, gj, di, gg, n)
    {
        init();
    }

    /**
     * @brief Узнать название предобуславливателя
     * @return Название предобуславливателя
     */
    virtual std::string get_name() const
    {
        return "Nothing_MKL";
    }

    /**
     * @brief Решает СЛАУ вида x = S^-1 * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_S(const val_type * f, val_type * x) const
    {
        copy_vec(f, x);
    }

    /**
     * @brief Решает СЛАУ вида x = S^-T * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_ST(const val_type * f, val_type * x) const
    {
        copy_vec(f, x);
    }

    /**
     * @brief Решает СЛАУ вида x = Q^-1 * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_Q(const val_type * f, val_type * x) const
    {
        copy_vec(f, x);
    }

    /**
     * @brief Решает СЛАУ вида x = Q^-T * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_QT(const val_type * f, val_type * x) const
    {
        copy_vec(f, x);
    }

    /**
     * @brief Вычисляет x = Q * f
     * @param[in] f
     * @param[out] x
     */
    virtual void mul_Q(const val_type * f, val_type * x) const
    {
        copy_vec(f, x);
    }

protected:

    void init()
    {
        wrappers::mkl::wrapper_mkl_set_env_max_threads();
    }

    inline void copy_vec(const float * x, float * y) const
    {
        cblas_scopy(static_cast<MKL_INT>(this->m_n), const_cast<float *>(x), 1, y, 1);
    }

    inline void copy_vec(const double * x, double * y) const
    {
        cblas_dcopy(static_cast<MKL_INT>(this->m_n), const_cast<double *>(x), 1, y, 1);
    }

    inline void copy_vec(const std::complex<float> * x, std::complex<float> * y) const
    {
        MKL_Complex8 * x1 = reinterpret_cast<MKL_Complex8 *>(const_cast<std::complex<float> *>(x));
        MKL_Complex8 * y1 = reinterpret_cast<MKL_Complex8 *>(const_cast<std::complex<float> *>(y));
        cblas_ccopy(static_cast<MKL_INT>(this->m_n), x1, 1, y1, 1);
    }

    void copy_vec(const std::complex<double> * x, std::complex<double> * y) const
    {
        MKL_Complex16 * x1 = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(x));
        MKL_Complex16 * y1 = reinterpret_cast<MKL_Complex16 *>(const_cast<std::complex<double> *>(y));
        cblas_zcopy(static_cast<MKL_INT>(this->m_n), x1, 1, y1, 1);
    }
};

}}}} // namespace core::solvers::CSRC::preconditioners

#endif // SOLVERS_CSRC_PRECONDITIONERS_NOTHING_MKL_H_INCLUDED
