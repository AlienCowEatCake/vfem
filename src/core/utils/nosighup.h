#if !defined(UTILS_NOSIGHUP_H_INCLUDED)
#define UTILS_NOSIGHUP_H_INCLUDED

namespace core { namespace utils { namespace nosighup {

/**
 * @brief Устанавливает обработчик SIGHUP с перенаправлением вывода в файлы
 * @note Hе использовать при отладке!
 */
void set_nosighup();

}}} // namespace core::utils::nosighup

#endif // UTILS_NOSIGHUP_H_INCLUDED
