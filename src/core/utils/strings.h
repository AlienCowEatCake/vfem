#if !defined(UTILS_STRINGS_H_INCLUDED)
#define UTILS_STRINGS_H_INCLUDED

#include <string>

namespace core { namespace utils { namespace strings {

/**
 * @brief Обрезка whitespace символов по краям строки
 * @param[in] str входная строка
 * @return обрезанная строка
 */
std::string trim(const std::string & str);

/**
 * @brief Переобразование к нижнему регистру (только для символов US-ASCII)
 * @param[in] str входная строка
 * @return строка в нижнем регистре
 */
std::string to_lowercase(const std::string & str);

/**
 * @brief Переобразование к верхнему регистру (только для символов US-ASCII)
 * @param[in] str входная строка
 * @return строка в верхнем регистре
 */
std::string to_uppercase(const std::string & str);

/**
 * @brief Выполняет сравнение строк без учета регистра (только для символов US-ASCII)
 * @param[in] str1 входная строка #1
 * @param[in] str2 входная строка #2
 * @return < 0 - str1 меньше str2, 0 - str1 идентично str2, > 0 - str1 больше str2
 */
int compare_ci(const std::string & str1, const std::string & str2);

}}} // namespace core::utils::strings

#endif // UTILS_STRINGS_H_INCLUDED
