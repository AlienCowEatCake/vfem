#include "solvers_factory.h"

#include "../../utils/strings.h"

#include "preconditioners/Nothing/preconditioner_Nothing.h"
#include "preconditioners/Nothing/preconditioner_Nothing_OpenMP.h"
#include "preconditioners/Nothing/preconditioner_Nothing_MKL.h"
#include "preconditioners/Di/preconditioner_Di.h"
#include "preconditioners/Di/preconditioner_Di_OpenMP.h"
#include "preconditioners/GS/preconditioner_GS.h"
#include "preconditioners/LDLT/preconditioner_LDLT.h"
#include "preconditioners/LLT/preconditioner_LLT.h"
#include "preconditioners/LLT/preconditioner_LLT_MKL.h"

#include "symmetric/complex/COCG/COCG.h"
#include "symmetric/complex/COCG/COCG_Smooth.h"
#include "symmetric/complex/COCG/COCG_Smooth_OpenMP.h"
#include "symmetric/complex/COCG/COCG_Smooth_MKL.h"

#include "symmetric/complex/COCR/COCR.h"
#include "symmetric/complex/COCR/COCR_Smooth.h"

#include "symmetric/complex/BiCG_Complex/BiCG_Complex.h"
#include "symmetric/complex/BiCG_Complex/BiCG_Complex_Smooth.h"

#include "symmetric/complex/BiCGStab_Complex/BiCGStab_Complex.h"
#include "symmetric/complex/BiCGStab_Complex/BiCGStab_Complex_Smooth.h"

#include "symmetric/complex/GMRES_Complex/GMRES_Complex.h"
#include "symmetric/complex/GMRES_Complex/GMRES_Complex_OpenMP.h"
#include "symmetric/complex/GMRES_Complex/GMRES_Complex_MKL.h"

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
        const std::complex<double> * di, const std::complex<double> * gg, std::size_t n)
{
    if(utils::strings::compare_ci(name, "Nothing") == 0)
        return new preconditioners::preconditioner_Nothing<std::complex<double>, std::size_t>(gi, gj, di, gg, n);
    if(utils::strings::compare_ci(name, "Nothing_OpenMP") == 0)
        return new preconditioners::preconditioner_Nothing_OpenMP<std::complex<double>, std::size_t>(gi, gj, di, gg, n);
    if(utils::strings::compare_ci(name, "Nothing_MKL") == 0)
        return new preconditioners::preconditioner_Nothing_MKL<std::complex<double>, std::size_t>(gi, gj, di, gg, n);
    if(utils::strings::compare_ci(name, "Di") == 0)
        return new preconditioners::preconditioner_Di<std::complex<double>, std::size_t>(gi, gj, di, gg, n);
    if(utils::strings::compare_ci(name, "Di_OpenMP") == 0)
        return new preconditioners::preconditioner_Di_OpenMP<std::complex<double>, std::size_t>(gi, gj, di, gg, n);
    if(utils::strings::compare_ci(name, "GS") == 0)
        return new preconditioners::preconditioner_GS<std::complex<double>, std::size_t>(gi, gj, di, gg, n);
    if(utils::strings::compare_ci(name, "LDLT") == 0)
        return new preconditioners::preconditioner_LDLT<std::complex<double>, std::size_t>(gi, gj, di, gg, n);
    if(utils::strings::compare_ci(name, "LLT") == 0)
        return new preconditioners::preconditioner_LLT<std::complex<double>, std::size_t>(gi, gj, di, gg, n);
    if(utils::strings::compare_ci(name, "LLT_MKL") == 0)
        return new preconditioners::preconditioner_LLT_MKL<std::size_t>(gi, gj, di, gg, n);
#if defined(USE_MKL)
    return new preconditioners::preconditioner_LLT_MKL<std::size_t>(gi, gj, di, gg, n);
#else
    return new preconditioners::preconditioner_LLT<std::complex<double>, std::size_t>(gi, gj, di, gg, n);
#endif
}

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
        preconditioners::preconditioner_interface<std::complex<double>, std::size_t> * precond)
{
    symmetric::symmetric_solver_interface<std::complex<double>, std::size_t> * solver;
    if(utils::strings::compare_ci(name, "COCG") == 0)
        solver = new symmetric::complex::COCG();
    else if(utils::strings::compare_ci(name, "COCG_Smooth") == 0)
        solver = new symmetric::complex::COCG_Smooth();
    else if(utils::strings::compare_ci(name, "COCG_Smooth_OpenMP") == 0)
        solver = new symmetric::complex::COCG_Smooth_OpenMP();
    else if(utils::strings::compare_ci(name, "COCG_Smooth_MKL") == 0)
        solver = new symmetric::complex::COCG_Smooth_MKL();
    else if(utils::strings::compare_ci(name, "COCR") == 0)
        solver = new symmetric::complex::COCR();
    else if(utils::strings::compare_ci(name, "COCR_Smooth") == 0)
        solver = new symmetric::complex::COCR_Smooth();
    else if(utils::strings::compare_ci(name, "BiCG_Complex") == 0 ||
            utils::strings::compare_ci(name, "BiCG") == 0)
        solver = new symmetric::complex::BiCG_Complex();
    else if(utils::strings::compare_ci(name, "BiCG_Complex_Smooth") == 0 ||
            utils::strings::compare_ci(name, "BiCG_Smooth") == 0)
        solver = new symmetric::complex::BiCG_Complex_Smooth();
    else if(utils::strings::compare_ci(name, "BiCGStab_Complex") == 0 ||
            utils::strings::compare_ci(name, "BiCGStab") == 0)
        solver = new symmetric::complex::BiCGStab_Complex();
    else if(utils::strings::compare_ci(name, "BiCGStab_Complex_Smooth") == 0 ||
            utils::strings::compare_ci(name, "BiCGStab_Smooth") == 0)
        solver = new symmetric::complex::BiCGStab_Complex_Smooth();
    else if(utils::strings::compare_ci(name, "GMRES_Complex") == 0 ||
            utils::strings::compare_ci(name, "GMRES") == 0)
        solver = new symmetric::complex::GMRES_Complex();
    else if(utils::strings::compare_ci(name, "GMRES_Complex_OpenMP") == 0 ||
            utils::strings::compare_ci(name, "GMRES_OpenMP") == 0)
        solver = new symmetric::complex::GMRES_Complex_OpenMP();
    else if(utils::strings::compare_ci(name, "GMRES_Complex_MKL") == 0 ||
            utils::strings::compare_ci(name, "GMRES_MKL") == 0)
        solver = new symmetric::complex::GMRES_Complex_MKL();
    else
#if defined(USE_MKL)
        solver = new symmetric::complex::COCG_Smooth_MKL();
#else
        solver = new symmetric::complex::COCG_Smooth();
#endif
    solver->init(gi, gj, di, gg, n, precond);
    return solver;
}

}}}} // namespace core::solvers::CSRC::factory
