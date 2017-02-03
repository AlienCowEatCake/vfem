#if !defined(SOLVERS_CSLR_PRECONDITIONERS_INTERFACE_H_INCLUDED)
#define SOLVERS_CSLR_PRECONDITIONERS_INTERFACE_H_INCLUDED

#include <cstddef>
#include <string>

namespace fem_core { namespace solvers { namespace CSLR { namespace preconditioners {

/**
 * @brief Интерфейс для предобуславливателей вида A = S * Q
 */
template<typename val_type, typename ind_type = std::size_t>
class preconditioner_interface
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
    preconditioner_interface(const ind_type * gi, const ind_type * gj, const val_type * di,
                             const val_type * gl, const val_type * gu, ind_type n);

    /**
     * @brief Конструктор для симметричного предобуславливателя
     * @param[in] gi Массив ig в разряженном строчно-столбцовом представлении
     * @param[in] gj Массив jg в разряженном строчно-столбцовом представлении
     * @param[in] di Массив диагональных эл-тов в разряженном строчно-столбцовом представлении
     * @param[in] gg Массив нижнего и верхнего треугольника в разряженном строчно-столбцовом представлении
     * @param[in] n Размерность матрицы
     */
    preconditioner_interface(const ind_type * gi, const ind_type * gj, const val_type * di,
                             const val_type * gg, ind_type n);

    /**
     * @brief Узнать название предобуславливателя
     * @return Название предобуславливателя
     */
    virtual std::string get_name() const = 0;

    /**
     * @brief Решает СЛАУ вида x = S^-1 * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_S(const val_type * f, val_type * x) const = 0;

    /**
     * @brief Решает СЛАУ вида x = S^-T * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_ST(const val_type * f, val_type * x) const = 0;

    /**
     * @brief Решает СЛАУ вида x = Q^-1 * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_Q(const val_type * f, val_type * x) const = 0;

    /**
     * @brief Решает СЛАУ вида x = Q^-T * f
     * @param[in] f
     * @param[out] x
     */
    virtual void solve_QT(const val_type * f, val_type * x) const = 0;

    /**
     * @brief Вычисляет x = Q * f
     * @param[in] f
     * @param[out] x
     */
    virtual void mul_Q(const val_type * f, val_type * x) const = 0;

    /**
     * @brief Виртуальный деструктор
     */
    virtual ~preconditioner_interface()
    {}

protected:
    preconditioner_interface()
    {}

private:
    preconditioner_interface(const preconditioner_interface & other);
    const preconditioner_interface & operator = (const preconditioner_interface & other);
};

}}}} // namespace fem_core::solvers::CSLR::preconditioners

#endif // SOLVERS_CSLR_PRECONDITIONERS_INTERFACE_H_INCLUDED
