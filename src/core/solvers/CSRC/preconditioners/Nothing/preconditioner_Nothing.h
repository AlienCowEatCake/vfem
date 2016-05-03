#if !defined(SOLVERS_CSRC_PRECONDITIONERS_NOTHING_H_INCLUDED)
#define SOLVERS_CSRC_PRECONDITIONERS_NOTHING_H_INCLUDED

#include "../preconditioner_interface.h"

namespace core { namespace solvers { namespace CSRC { namespace preconditioners {

/**
 * @brief Класс базовый (пустой) предобуславливатель A = S * Q
 */
template<typename val_type, typename ind_type = std::size_t>
class preconditioner_Nothing : public preconditioner_interface<val_type, ind_type>
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
    preconditioner_Nothing(const ind_type * gi, const ind_type * gj, const val_type * di,
                           const val_type * gl, const val_type * gu, ind_type n)
        : m_gi(gi), m_gj(gj), m_di(di), m_gl(gl), m_gu(gu), m_n(n), m_is_symmetric(false)
    {}

    /**
     * @brief Конструктор для симметричного предобуславливателя
     * @param[in] gi Массив ig в разряженном строчно-столбцовом представлении
     * @param[in] gj Массив jg в разряженном строчно-столбцовом представлении
     * @param[in] di Массив диагональных эл-тов в разряженном строчно-столбцовом представлении
     * @param[in] gg Массив нижнего и верхнего треугольника в разряженном строчно-столбцовом представлении
     * @param[in] n Размерность матрицы
     */
    preconditioner_Nothing(const ind_type * gi, const ind_type * gj, const val_type * di,
                           const val_type * gg, ind_type n)
        : m_gi(gi), m_gj(gj), m_di(di), m_gl(gg), m_gu(gg), m_n(n), m_is_symmetric(true)
    {}

    /**
     * @brief Узнать название предобуславливателя
     * @return Название предобуславливателя
     */
    virtual std::string get_name() const
    {
        return "Nothing";
    }

    /**
     * @brief Решает СЛАУ вида x = S^-1 * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_S(const val_type * f, val_type * x) const
    {
        for(ind_type k = 0; k < m_n; k++)
            x[k] = f[k];
    }

    /**
     * @brief Решает СЛАУ вида x = S^-T * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_ST(const val_type * f, val_type * x) const
    {
        for(ind_type k = 0; k < m_n; k++)
            x[k] = f[k];
    }

    /**
     * @brief Решает СЛАУ вида x = Q^-1 * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_Q(const val_type * f, val_type * x) const
    {
        for(ind_type k = 0; k < m_n; k++)
            x[k] = f[k];
    }

    /**
     * @brief Решает СЛАУ вида x = Q^-T * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_QT(const val_type * f, val_type * x) const
    {
        for(ind_type k = 0; k < m_n; k++)
            x[k] = f[k];
    }

    /**
     * @brief Вычисляет x = Q * f
     * @param[in] f
     * @param[out] x
     */
    virtual void mul_Q(const val_type * f, val_type * x) const
    {
        for(ind_type k = 0; k < m_n; k++)
            x[k] = f[k];
    }

    /**
     * @brief Виртуальный деструктор
     */
    virtual ~preconditioner_Nothing()
    {}

protected:

    const ind_type * m_gi, * m_gj;
    const val_type * m_di, * m_gl, * m_gu;
    ind_type m_n;
    bool m_is_symmetric;
};

}}}} // namespace core::solvers::CSRC::preconditioners

#endif // SOLVERS_CSRC_PRECONDITIONERS_NOTHING_H_INCLUDED
