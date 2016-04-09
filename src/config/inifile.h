#ifndef INIFILE_H
#define INIFILE_H

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <list>
#include <set>

// Обрезка whitespace символов по краям строки
std::string trim(const std::string & str);

// Переобразование к нижнему регистру (только для символов US-ASCII)
std::string to_lowercase(const std::string & str);

class inifile
{
public:
    // Немного конструкторов
    inifile() : status(false) {}
    inifile(const std::string & filename) { status = load(filename); }
    // Загрузка ini-файла с именем "filename"
    bool load(const std::string & filename);
    // Статус, разобран ли ini-файл
    inline bool good() const { return status; }

    // Функция для типа string, string не стоит получать с помощью stringstream, так как может многое потеряться
    // Получить значение параметра с именем "parameter" в секции "section" подсекции "subsection".
    // Если такой параметр не найден, вернуть значение по-умолчанию "fallback".
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

    // Функция для типа bool, у него нужно сделать дополнительное преобразование yes|true|1 -> true
    // Получить значение параметра с именем "parameter" в секции "section" подсекции "subsection".
    // Если такой параметр не найден, вернуть значение по-умолчанию "fallback".
    template<typename U>
    bool get(const std::string & section, const U & subsection, const std::string & parameter, const bool & fallback) const
    {
        // Попробуем найти искомое значение
        const std::string * result_pstr = get_internal(section, subsection_to_str(subsection), parameter);
        if(result_pstr)
        {
            std::string value = to_lowercase(* result_pstr);
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

    // Обобщенная функция для типов, которые не являются типами bool или string
    // Получить значение параметра с именем "parameter" в секции "section" подсекции "subsection".
    // Если такой параметр не найден, вернуть значение по-умолчанию "fallback".
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

    // Получть список из всех подсекций секции "section". Тип возвращаемого результата определяется
    // типом указателя "type". Сам указатель является фиктивным и нигде не используется, то есть
    // можно туда передавать что угодно, хоть (size_t)(NULL), например.
    template<typename T>
    std::list<T> enumerate(const std::string & section, const T * type) const
    {
        (void)(type);
        std::list<T> result;
        typename std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::const_iterator it_section = values.find(to_lowercase(section));
        if(it_section != values.end())
        {
            for(std::map<std::string, std::map<std::string, std::string> >::const_iterator it = it_section->second.begin(); it != it_section->second.end(); ++it)
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

    // Проверить, что в файле нет никаких секций, кроме секций из списка "whitelist"
    // Полезно для контроля опечаток, так как в остальных местах это штатная ситуация
    bool check_sections(const std::list<std::string> & whitelist) const;

    // Проверить, что в секции "section" во всех подсекциях нет никаких параметров, кроме параметров из списка "whitelist"
    // Полезно для контроля опечаток, так как в остальных местах это штатная ситуация
    bool check_parameters(const std::string & section, const std::list<std::string> & whitelist) const;

protected:

    // Статус, разобран ли входной конфиг
    bool status;

    // Конетейнер значений параметров, сгруппированных в подсекции, сгруппированные в секции
    std::map<std::string, std::map<std::string, std::map<std::string, std::string> > > values;

    // Получить значение параметра с именем "parameter" в секции "section" подсекции "subsection".
    const std::string * get_internal(const std::string & section, const std::string & subsection, const std::string & parameter) const;

    // Преобразование подсекции в строковое представление
    // Если это строка, то она уже итак задана как надо
    inline std::string subsection_to_str(const std::string & subsection) const
    {
        return subsection;
    }

    // Преобразование подсекции в строковое представление
    // Если это char *, то просто соберем из него строку
    inline std::string subsection_to_str(const char * subsection) const
    {
        return std::string(subsection);
    }

    // Преобразование подсекции в строковое представление
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


#endif // INIFILE_H
