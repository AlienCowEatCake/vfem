#include "omp_stubs.h"

namespace omp_stubs
{

/**
 * @brief Sets the number of threads in subsequent parallel regions, unless overridden by a num_threads clause.
 * @param[in] The number of threads in the parallel region.
 */
void omp_set_num_threads(int num_threads)
{
    (void)num_threads;
}

/**
 * @brief Returns the number of threads in the parallel region.
 * @return The number of threads in the parallel region.
 */
int omp_get_num_threads()
{
    return 1;
}

/**
 * @brief Returns the thread number of the thread executing within its thread team.
 * @return The thread number of the thread executing within its thread team.
 */
int omp_get_thread_num()
{
    return 1;
}

/**
 * @brief Returns an integer that is equal to or greater than the number of threads that would be available if a parallel region without num_threads were defined at that point in the code.
 * @return An integer that is equal to or greater than the number of threads that would be available if a parallel region without num_threads were defined at that point in the code.
 */
int omp_get_max_threads()
{
    return 1;
}

}
