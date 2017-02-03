#include "strings.h"
#include <algorithm>
#include <cctype>

namespace fem_core { namespace utils { namespace strings {

/**
 * @brief Обрезка whitespace символов по краям строки
 * @param[in] str входная строка
 * @return обрезанная строка
 */
std::string trim(const std::string & str)
{
    size_t start = str.find_first_not_of(" \t\f\v\n\r");
    size_t stop = str.find_last_not_of(" \t\f\v\n\r");
    if(start == std::string::npos || stop == std::string::npos)
        return "";
    return str.substr(start, stop - start + 1);
}

/**
 * @brief Переобразование к нижнему регистру (только для символов US-ASCII)
 * @param[in] str входная строка
 * @return строка в нижнем регистре
 */
std::string to_lowercase(const std::string & str)
{
    std::string result = str;
    std::transform(str.begin(), str.end(), result.begin(), ::tolower);
    return result;
}

/**
 * @brief Переобразование к верхнему регистру (только для символов US-ASCII)
 * @param[in] str входная строка
 * @return строка в верхнем регистре
 */
std::string to_uppercase(const std::string & str)
{
    std::string result = str;
    std::transform(str.begin(), str.end(), result.begin(), ::toupper);
    return result;
}

/**
 * @brief Выполняет сравнение строк без учета регистра (только для символов US-ASCII)
 * @param[in] str1 входная строка #1
 * @param[in] str2 входная строка #2
 * @return < 0 - str1 меньше str2, 0 - str1 идентично str2, > 0 - str1 больше str2
 */
int compare_ci(const std::string & str1, const std::string & str2)
{
    return to_lowercase(str1).compare(to_lowercase(str2));
}

}}} // namespace fem_core::utils::strings
