#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "omp_wrapper.h"
#include <cstdlib>

// *************************************************************************************************

#if !defined(USE_OMP)

namespace fem_core { namespace wrappers { namespace omp_stubs {

/**
 * @brief Устанавливает количество потоков в последующих параллельных регионах.
 * @param[in] Количество потоков в параллельном регионе.
 */
void omp_set_num_threads(int num_threads)
{
    (void)num_threads;
}

/**
 * @brief Возвращает количество потоков в параллельном регионе.
 * @return Количество потоков в параллельном регионе.
 */
int omp_get_num_threads()
{
    return 1;
}

/**
 * @brief Возвращает номер выполняемого потока.
 * @return Номер выполняемого потока.
 */
int omp_get_thread_num()
{
    return 1;
}

/**
 * @brief Возвращает целое число, которое больше или равно количеству потоков для параллельных регионов.
 * @return Целое число, которое больше или равно количеству потоков для параллельных регионов.
 */
int omp_get_max_threads()
{
    return 1;
}

}}} // namespace fem_core::wrappers::omp_stubs

#endif

// *************************************************************************************************

namespace fem_core { namespace wrappers { namespace omp {

/**
 * @brief Установить количество тредов для OpenMP равным переменной окружения OMP_NUM_THREADS
 */
void wrapper_omp_set_env_max_threads()
{
    int num_threads = -1;
    if(num_threads <= 0)
    {
        char * env_num_threads = getenv("OMP_NUM_THREADS");
        if(env_num_threads)
            num_threads = atoi(env_num_threads);
    }
    if(num_threads <= 0)
    {
        num_threads = 1;
    }
    omp_set_num_threads(num_threads);
}

}}} // namespace fem_core::wrappers::omp

// *************************************************************************************************
