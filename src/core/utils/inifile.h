#if !defined(UTILS_INIFILE_H_INCLUDED)
#define UTILS_INIFILE_H_INCLUDED

#include <list>
#include <map>
#include <string>
#include <sstream>
#include <iostream>

#include "strings.h"

namespace core { namespace utils {

/**
 * @brief Класс для работы с ini-файлами
 */
class inifile
{
public:
    /**
     * @brief Конструктор по-умолчанию
     */
    inifile() : m_status(false) {}

    /**
     * @brief Конструктор с именем входного файла
     * @param[in] filename имя входного файла
     */
    inifile(const std::string & filename) { m_status = load(filename); }

    /**
     * @brief Загрузка ini-файла
     * @param[in] filename имя входного файла
     * @return true - файл успешно считан, false - возникли ошибки
     */
    bool load(const std::string & filename);

    /**
     * @brief Статус, разобран ли ini-файл
     * @return true - файл успешно разобран, false - возникли ошибки
     */
    inline bool good() const { return m_status; }

    /**
     * @brief Получить значение параметра с именем "parameter" в секции "section" подсекции "subsection"
     * @param[in] section секция, в которой нужно искать
     * @param[in] subsection подсекция, в которой нужно искать
     * @param[in] parameter параметр, который нужно искать
     * @param[in] fallback значение по-умолчанию, которое возвращается, если искомый параметр не найден
     * @return значение параметра или fallback, если параметр не найден
     * @note Функция для типа string, string не стоит получать с помощью stringstream, так как может многое потеряться
     */
    template<typename U>
    std::string get(const std::string & section, const U & subsection, const std::string & parameter, const std::string & fallback) const
    {
        // Попробуем найти искомое значение
        const std::string * result_pstr = get_internal(section, subsection_to_str(subsection), parameter);
        if(result_pstr)
        {
            std::string result = * result_pstr;
            //std::cout << section << "." << subsection_to_str(subsection) << " : " << parameter << " = " << result << std::endl;
            return result;
        }
        // Ничего не нашли, вернем значение по-умолчанию
        return fallback;
    }

    /**
     * @brief Получить значение параметра с именем "parameter" в секции "section" подсекции "subsection"
     * @param[in] section секция, в которой нужно искать
     * @param[in] subsection подсекция, в которой нужно искать
     * @param[in] parameter параметр, который нужно искать
     * @param[in] fallback значение по-умолчанию, которое возвращается, если искомый параметр не найден
     * @return значение параметра или fallback, если параметр не найден
     * @note Функция для типа char *, его приводим к std::string
     */
    template<typename U>
    inline std::string get(const std::string & section, const U & subsection, const std::string & parameter, const char * fallback) const
    {
        return get(section, subsection, parameter, std::string(fallback));
    }

    /**
     * @brief Получить значение параметра с именем "parameter" в секции "section" подсекции "subsection"
     * @param[in] section секция, в которой нужно искать
     * @param[in] subsection подсекция, в которой нужно искать
     * @param[in] parameter параметр, который нужно искать
     * @param[in] fallback значение по-умолчанию, которое возвращается, если искомый параметр не найден
     * @return значение параметра или fallback, если параметр не найден
     * @note Функция для типа bool, у него нужно сделать дополнительное преобразование yes|true|1 -> true
     */
    template<typename U>
    bool get(const std::string & section, const U & subsection, const std::string & parameter, const bool & fallback) const
    {
        // Попробуем найти искомое значение
        const std::string * result_pstr = get_internal(section, subsection_to_str(subsection), parameter);
        if(result_pstr)
        {
            std::string value = strings::to_lowercase(* result_pstr);
            bool result;
            if(value == "yes" || value == "true" || value == "1")
                result = true;
            else
                result = false;
            //std::cout << section << "." << subsection_to_str(subsection) << " : " << parameter << " = " << result << std::endl;
            return result;
        }
        // Ничего не нашли, вернем значение по-умолчанию
        return fallback;
    }

    /**
     * @brief Получить значение параметра с именем "parameter" в секции "section" подсекции "subsection"
     * @param[in] section секция, в которой нужно искать
     * @param[in] subsection подсекция, в которой нужно искать
     * @param[in] parameter параметр, который нужно искать
     * @param[in] fallback значение по-умолчанию, которое возвращается, если искомый параметр не найден
     * @return значение параметра или fallback, если параметр не найден
     * @note Обобщенная функция для типов, которые не являются типами bool или string
     */
    template<typename T, typename U>
    T get(const std::string & section, const U & subsection, const std::string & parameter, const T & fallback) const
    {
        // Попробуем найти искомое значение
        const std::string * result_pstr = get_internal(section, subsection_to_str(subsection), parameter);
        if(result_pstr)
        {
            std::stringstream sst(* result_pstr);
            T result;
            sst >> result;
            //std::cout << section << "." << subsection_to_str(subsection) << " : " << parameter << " = " << result << std::endl;
            return result;
        }
        // Ничего не нашли, вернем значение по-умолчанию
        return fallback;
    }

    /**
     * @brief Получть список из всех подсекций секции "section"
     * @param[in] section секция, в которой нужно искать подсекции
     * @return список из всех подсекций секции "section"
     */
    template<typename T>
    std::list<T> enumerate(const std::string & section) const
    {
        std::list<T> result;
        std::map<std::string, section_t>::const_iterator it_section = m_values.find(strings::to_lowercase(section));
        if(it_section != m_values.end())
        {
            for(std::map<std::string, subsection_t>::const_iterator
                it = it_section->second._.begin(), it_end = it_section->second._.end(); it != it_end; ++it)
            {
                std::stringstream sst(it->first);
                T subsection = T();
                sst >> subsection;
                //std::cout << " * enumerate(" << section << "): " << subsection << std::endl;
                result.push_back(subsection);
            }
        }
        return result;
    }

    /**
     * @brief Проверить, что в файле нет никаких секций, кроме секций из списка "whitelist"
     * @param[in] whitelist список разрешенных секций
     * @return true - все секции присутствуют в "whitelist", false - есть неразрешенная секция
     * @note Полезно для контроля опечаток, так как в остальных местах это штатная ситуация
     */
    bool check_sections(const std::list<std::string> & whitelist) const;

    /**
     * @brief Проверить, что в секции "section" во всех подсекциях нет никаких параметров, кроме параметров из списка "whitelist"
     * @param[in] section секция, для которой будет выполнена проверка
     * @param[in] whitelist список разрешенных параметров
     * @return true - все параметры присутствуют в "whitelist", false - есть неразрешенный параметр
     * @note Полезно для контроля опечаток, так как в остальных местах это штатная ситуация
     */
    bool check_parameters(const std::string & section, const std::list<std::string> & whitelist) const;

protected:

    /**
     * @brief Костыльный контейнер подсекции для исправления багов MSVC (https://msdn.microsoft.com/en-us/library/074af4b6.aspx)
     */
    struct subsection_t
    {
        std::map<std::string, std::string> _;
    };

    /**
     * @brief Костыльный контейнер секции для исправления багов MSVC (https://msdn.microsoft.com/en-us/library/074af4b6.aspx)
     */
    struct section_t
    {
        std::map<std::string, subsection_t> _;
    };

    /**
     * @brief Статус, разобран ли входной конфиг
     */
    bool m_status;

    /**
     * @brief Конетейнер значений параметров, сгруппированных в подсекции, сгруппированные в секции
     */
    std::map<std::string, section_t> m_values;

    /**
     * @brief Получить значение параметра с именем "parameter" в секции "section" подсекции "subsection"
     * @param[in] section секция
     * @param[in] subsection подсекция
     * @param[in] parameter параметр
     * @return значение параметра или NULL, если такого параметра нет
     */
    const std::string * get_internal(const std::string & section, const std::string & subsection, const std::string & parameter) const;

    /**
     * @brief Преобразование подсекции в строковое представление
     * @param[in] subsection подсекция
     * @return строковое представление подсекции
     * @note Если это строка, то она уже итак задана как надо
     */
    inline std::string subsection_to_str(const std::string & subsection) const
    {
        return subsection;
    }

    /**
     * @brief Преобразование подсекции в строковое представление
     * @param[in] subsection подсекция
     * @return строковое представление подсекции
     * @note Если это char *, то просто соберем из него строку
     */
    inline std::string subsection_to_str(const char * subsection) const
    {
        return std::string(subsection);
    }

    /**
     * @brief Преобразование подсекции в строковое представление
     * @param[in] subsection подсекция
     * @return строковое представление подсекции
     */
    template<typename T>
    std::string subsection_to_str(const T & subsection) const
    {
        std::string result;
        // Проверим, эквивалентна ли эта строка пустой подсекции, это нужно,
        // например, для int-типов, так как для них пустая подсекция - "0"
        std::stringstream sst("");
        T subsection_test = T();
        sst >> subsection_test;
        // Если не эквивалентна, то нужно выполнить преобразование
        if(subsection_test != subsection)
        {
            // Преобразуем подсекцию в строку
            std::stringstream subsection_sst;
            subsection_sst << subsection;
            result = subsection_sst.str();
        }
        //std::cout << " * convert(\"" << subsection << "\") -> \"" << result << "\"" << std::endl;
        return result;
    }
};

}} // namespace core::utils

#endif // UTILS_INIFILE_H_INCLUDED
