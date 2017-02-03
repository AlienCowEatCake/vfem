#if !defined(UTILS_TIMERS_H_INCLUDED)
#define UTILS_TIMERS_H_INCLUDED

#include <string>

namespace fem_core { namespace utils { namespace timers {

/**
 * @brief Возвращает время в миллисекундах, само по себе ничего не значит, но по разнице можно засекать время
 * @return Время в миллисекундах
 */
unsigned long mtime();

/**
 * @brief Функция для распечатки затраченного времени, предназначена для работы совместно с mtime()
 * @param[in] msec Затраченное время в миллисекундах
 * @param[in] descr Текстовое описание операции
 */
void print_time(unsigned long msec, const std::string & descr);

}}} // namespace fem_core::utils::timers

#endif // UTILS_TIMERS_H_INCLUDED
