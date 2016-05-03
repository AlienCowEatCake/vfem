#if !defined(SOLVERS_SYMMETRIC_INTERFACE_H_INCLUDED)
#define SOLVERS_SYMMETRIC_INTERFACE_H_INCLUDED

#include <cstdlib>

#include "../preconditioners/preconditioner_interface.h"

namespace core { namespace solvers { namespace CSRC { namespace symmetric {

/**
 * @brief Интерфейс для симметричных решателей
 */
template<typename val_type, typename ind_type = std::size_t>
class symmetric_solver_interface
{
public:

    /**
     * @brief Инициализация симметричная
     * @param[in] gi Массив ig в разряженном строчно-столбцовом представлении
     * @param[in] gj Массив jg в разряженном строчно-столбцовом представлении
     * @param[in] di Массив диагональных эл-тов в разряженном строчно-столбцовом представлении
     * @param[in] gg Массив нижнего и верхнего треугольника в разряженном строчно-столбцовом представлении
     * @param[in] n Размерность матрицы
     * @param[in] precond Предобуславливатель
     */
    virtual void init(const ind_type * gi, const ind_type * gj, const val_type * di,
                      const val_type * gg, ind_type n,
                      preconditioners::preconditioner_interface<val_type, ind_type> * precond) = 0;

    /**
     * @brief Запуск решения
     * @param[inuot] solution На вход идет начальное приближение, на выход - решение
     * @param[in] rp Правая часть
     * @param[in] eps Целевая относительная невязка
     * @param[in] max_iter Максимальное количество итераций
     */
    virtual void solve(val_type * solution, const val_type * rp,
                       double eps, ind_type max_iter) = 0;

    /**
     * @brief Виртуальный деструктор
     */
    virtual ~symmetric_solver_interface() {}
};

}}}} // namespace core::solvers::CSRC::symmetric

#endif // SOLVERS_SYMMETRIC_INTERFACE_H_INCLUDED

