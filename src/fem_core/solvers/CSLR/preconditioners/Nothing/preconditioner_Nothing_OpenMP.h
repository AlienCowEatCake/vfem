#if !defined(SOLVERS_CSLR_PRECONDITIONERS_NOTHING_OPENMP_H_INCLUDED)
#define SOLVERS_CSLR_PRECONDITIONERS_NOTHING_OPENMP_H_INCLUDED

#include "preconditioner_Nothing.h"
#include "../../../../wrappers/omp_wrapper.h"

namespace fem_core { namespace solvers { namespace CSLR { namespace preconditioners {

/**
 * @brief Класс базовый (пустой) предобуславливатель A = S * Q, распараллеленный через OpenMP
 */
template<typename val_type, typename ind_type = std::size_t>
class preconditioner_Nothing_OpenMP : public preconditioner_Nothing<val_type, ind_type>
{
    typedef wrappers::omp::omp_int omp_int;

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
    preconditioner_Nothing_OpenMP(const ind_type * gi, const ind_type * gj, const val_type * di,
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
    preconditioner_Nothing_OpenMP(const ind_type * gi, const ind_type * gj, const val_type * di,
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
        return "Nothing_OpenMP";
    }

    /**
     * @brief Решает СЛАУ вида x = S^-1 * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_S(const val_type * f, val_type * x) const
    {
#pragma omp parallel for
        for(omp_int k = 0; k < static_cast<omp_int>(this->m_n); k++)
            x[k] = f[k];
    }

    /**
     * @brief Решает СЛАУ вида x = S^-T * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_ST(const val_type * f, val_type * x) const
    {
#pragma omp parallel for
        for(omp_int k = 0; k < static_cast<omp_int>(this->m_n); k++)
            x[k] = f[k];
    }

    /**
     * @brief Решает СЛАУ вида x = Q^-1 * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_Q(const val_type * f, val_type * x) const
    {
#pragma omp parallel for
        for(omp_int k = 0; k < static_cast<omp_int>(this->m_n); k++)
            x[k] = f[k];
    }

    /**
     * @brief Решает СЛАУ вида x = Q^-T * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_QT(const val_type * f, val_type * x) const
    {
#pragma omp parallel for
        for(omp_int k = 0; k < static_cast<omp_int>(this->m_n); k++)
            x[k] = f[k];
    }

    /**
     * @brief Вычисляет x = Q * f
     * @param[in] f
     * @param[out] x
     */
    virtual void mul_Q(const val_type * f, val_type * x) const
    {
#pragma omp parallel for
        for(omp_int k = 0; k < static_cast<omp_int>(this->m_n); k++)
            x[k] = f[k];
    }

protected:

    void init()
    {
        wrappers::omp::wrapper_omp_set_env_max_threads();
    }
};

}}}} // namespace fem_core::solvers::CSLR::preconditioners

#endif // SOLVERS_CSLR_PRECONDITIONERS_NOTHING_OPENMP_H_INCLUDED
