#if !defined(UTILS_PROGRESS_H_INCLUDED)
#define UTILS_PROGRESS_H_INCLUDED

#include <cstddef>

namespace fem_core { namespace utils { namespace progress {

/**
 * @brief Распечатывает на stdout стильный прогрессбар
 * @param[in] message Заголовок прогрессбара (может быть пустым)
 * @param[in] index Индекс текущего обрабатываемого элемента
 * @param[in] num_all Общее количество обрабатываемых элементов
 */
void show_progress(const char * message, std::size_t index, std::size_t num_all);

}}} // namespace fem_core::utils::progress

#endif // UTILS_PROGRESS_H_INCLUDED
