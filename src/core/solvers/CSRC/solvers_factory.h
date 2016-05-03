#if !defined(SOLVERS_CSRC_FACTORY_H_INCLUDED)
#define SOLVERS_CSRC_FACTORY_H_INCLUDED

#include <complex>
#include <string>

#include "preconditioners/preconditioner_interface.h"
#include "symmetric/symmetric_solver_interface.h"

namespace core { namespace solvers { namespace CSRC { namespace factory {

/**
 * @brief Функция, создающая новый симметричный комплексный предобуславливатель
 * @param[in] name Название предобуславливателя
 * @param[in] gi Массив ig в разряженном строчно-столбцовом представлении
 * @param[in] gj Массив jg в разряженном строчно-столбцовом представлении
 * @param[in] di Массив диагональных эл-тов в разряженном строчно-столбцовом представлении
 * @param[in] gg Массив нижнего и верхнего треугольника в разряженном строчно-столбцовом представлении
 * @param[in] n Размерность матрицы
 * @return Новый симметричный комплексный предобуславливатель (удаляет вызывающая сторона)
 */
preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * create_symmetric_complex_preconditioner(
        const std::string & name, const std::size_t * gi, const std::size_t * gj,
        const std::complex<double> * di, const std::complex<double> * gg, std::size_t n);

/**
 * @brief Функция, создающая новый симметричный комплексный решатель
 * @param[in] name Название решателя
 * @param[in] gi Массив ig в разряженном строчно-столбцовом представлении
 * @param[in] gj Массив jg в разряженном строчно-столбцовом представлении
 * @param[in] di Массив диагональных эл-тов в разряженном строчно-столбцовом представлении
 * @param[in] gg Массив нижнего и верхнего треугольника в разряженном строчно-столбцовом представлении
 * @param[in] n Размерность матрицы
 * @param[in] precond Предобуславливатель
 * @return Новый симметричный комплексный решатель (удаляет вызывающая сторона)
 */
symmetric::symmetric_solver_interface<std::complex<double>, std::size_t> * create_symmetric_complex_solver(
        const std::string & name, const std::size_t * gi, const std::size_t * gj,
        const std::complex<double> * di, const std::complex<double> * gg, std::size_t n,
        preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond);

}}}} // namespace core::solvers::CSRC::factory

#endif // SOLVERS_CSRC_FACTORY_H_INCLUDED
