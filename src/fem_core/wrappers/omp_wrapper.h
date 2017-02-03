#if !defined(WRAPPERS_OMP_WRAPPER_H_INCLUDED)
#define WRAPPERS_OMP_WRAPPER_H_INCLUDED

// *************************************************************************************************

#if defined(USE_OMP)

#include <omp.h>

#else

namespace fem_core { namespace wrappers { namespace omp_stubs {

/**
 * @brief Устанавливает количество потоков в последующих параллельных регионах.
 * @param[in] Количество потоков в параллельном регионе.
 */
void omp_set_num_threads(int num_threads);

/**
 * @brief Возвращает количество потоков в параллельном регионе.
 * @return Количество потоков в параллельном регионе.
 */
int omp_get_num_threads();

/**
 * @brief Возвращает номер выполняемого потока.
 * @return Номер выполняемого потока.
 */
int omp_get_thread_num();

/**
 * @brief Возвращает целое число, которое больше или равно количеству потоков для параллельных регионов.
 * @return Целое число, которое больше или равно количеству потоков для параллельных регионов.
 */
int omp_get_max_threads();

}}} // namespace fem_core::wrappers::omp_stubs

using namespace fem_core::wrappers::omp_stubs;

#endif

// *************************************************************************************************

#include <cstddef>

namespace fem_core { namespace wrappers { namespace omp {

#if defined(_MSC_VER)
typedef long omp_int;
#else
typedef std::size_t omp_int;
#endif

/**
 * @brief Установить количество тредов для OpenMP равным переменной окружения OMP_NUM_THREADS
 */
void wrapper_omp_set_env_max_threads();

}}} // namespace fem_core::wrappers::omp

// *************************************************************************************************

#endif // WRAPPERS_OMP_WRAPPER_H_INCLUDED
