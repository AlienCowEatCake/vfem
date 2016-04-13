#if !defined(OMP_STUBS_H)
#define OMP_STUBS_H

namespace omp_stubs
{

/**
 * @brief Sets the number of threads in subsequent parallel regions, unless overridden by a num_threads clause.
 * @param[in] The number of threads in the parallel region.
 */
void omp_set_num_threads(int num_threads);

/**
 * @brief Returns the number of threads in the parallel region.
 * @return The number of threads in the parallel region.
 */
int omp_get_num_threads();

/**
 * @brief Returns the thread number of the thread executing within its thread team.
 * @return The thread number of the thread executing within its thread team.
 */
int omp_get_thread_num();

}

#endif // OMP_STUBS_H
