#include "inifile.h"
#include <fstream>
#include <algorithm>
#include <cctype>
#include <set>

/**
 * @brief Обрезка whitespace символов по краям строки
 * @param str входная строка
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
 * @param str входная строка
 * @return строка в нижнем регистре
 */
std::string to_lowercase(const std::string & str)
{
    std::string result = str;
    std::transform(str.begin(), str.end(), result.begin(), ::tolower);
    return result;
}


/**
 * @brief Загрузка ini-файла
 * @param filename имя входного файла
 * @return true - файл успешно считан, false - возникли ошибки
 */
bool inifile::load(const std::string & filename)
{
    std::ifstream ifs(filename.c_str());
    if(!ifs.good())
    {
        status = false;
        return false;
    }

    std::string line;
    std::getline(ifs, line);
    line = trim(line);
    do
    {
        if(line.length() > 1 && line[0] == '[')
        {
            line = to_lowercase(line.substr(1, line.length() - 2));
            std::string section = line, subsection;
            size_t dot_pos = line.find_first_of(".");
            if(dot_pos != std::string::npos)
            {
                section = line.substr(0, dot_pos);
                subsection = line.substr(dot_pos + 1);
            }

            do
            {
                getline(ifs, line);
                line = trim(line);
                if(line.length() > 1 && line[0] != ';')
                {
                    size_t eq_pos = line.find_first_of("=");
                    if(eq_pos != std::string::npos)
                    {
                        std::string param = to_lowercase(trim(line.substr(0, eq_pos)));
                        std::string value = trim(line.substr(eq_pos + 1));
                        if(value.length() > 1 && value[0] == '\"')
                            value = trim(value.substr(1, value.length() - 2));
                        values[section][subsection][param] = value;
                    }
                }
            }
            while(ifs.good() && !(line.length() > 1 && line[0] == '['));
        }
        else
        {
            getline(ifs, line);
            line = trim(line);
        }
    }
    while(ifs.good());

    ifs.close();
    status = true;
    return true;
}

/**
 * @brief Получить значение параметра с именем "parameter" в секции "section" подсекции "subsection"
 * @param section секция
 * @param subsection подсекция
 * @param parameter параметр
 * @return значение параметра или NULL, если такого параметра нет
 */
const std::string * inifile::get_internal(const std::string & section, const std::string & subsection, const std::string & parameter) const
{
    // Найдем искомую секцию
    std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::const_iterator it_section = values.find(to_lowercase(section));
    if(it_section != values.end())
    {
        // Найдем искомую подсекцию
        std::map<std::string, std::map<std::string, std::string> >::const_iterator it_subsection = it_section->second.find(to_lowercase(subsection));
        if(it_subsection != it_section->second.end())
        {
            // И найдем искомый параметр
            std::map<std::string, std::string>::const_iterator it_parameter = it_subsection->second.find(to_lowercase(parameter));
            if(it_parameter != it_subsection->second.end())
            {
                return & (it_parameter->second);
            }
        }
    }
    // Ничего не нашли
    return NULL;
}

namespace
{

/**
 * @brief Преобразование std::list<std::string> -> std::set<std::string> c приведением к нижнему регистру
 * @param in входной std::list<std::string>
 * @return выходной std::set<std::string>
 */
std::set<std::string> list_to_set_and_lower(const std::list<std::string> & in)
{
    std::set<std::string> out;
    for(std::list<std::string>::const_iterator it = in.begin(), it_end = in.end(); it != it_end; ++it)
        out.insert(to_lowercase(* it));
    return out;
}

}

/**
 * @brief Проверить, что в файле нет никаких секций, кроме секций из списка "whitelist"
 * @param whitelist список разрешенных секций
 * @return true - все секции присутствуют в "whitelist", false - есть неразрешенная секция
 * @note Полезно для контроля опечаток, так как в остальных местах это штатная ситуация
 */
bool inifile::check_sections(const std::list<std::string> & whitelist) const
{
    if(!status)
        return false;
    std::set<std::string> whiteset = list_to_set_and_lower(whitelist);
    for(std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::const_iterator
        it = values.begin(), it_end = values.end(); it != it_end; ++it)
    {
        if(whiteset.find(it->first) == whiteset.end())
        {
            std::cerr << "[IniFile] Not-whitelisted section \"" << it->first << "\"" << std::endl;
            return false;
        }
    }
    return true;
}

/**
 * @brief Проверить, что в секции "section" во всех подсекциях нет никаких параметров, кроме параметров из списка "whitelist"
 * @param section секция, для которой будет выполнена проверка
 * @param whitelist список разрешенных параметров
 * @return true - все параметры присутствуют в "whitelist", false - есть неразрешенный параметр
 * @note Полезно для контроля опечаток, так как в остальных местах это штатная ситуация
 */
bool inifile::check_parameters(const std::string & section, const std::list<std::string> & whitelist) const
{
    if(!status)
        return false;
    std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::const_iterator
            it_section = values.find(to_lowercase(section));
    if(it_section == values.end())
        return true;
    std::set<std::string> whiteset = list_to_set_and_lower(whitelist);
    for(std::map<std::string, std::map<std::string, std::string> >::const_iterator it_subsection = it_section->second.begin(),
        it_subsection_end = it_section->second.end(); it_subsection != it_subsection_end; ++it_subsection)
    {
        for(std::map<std::string, std::string>::const_iterator it = it_subsection->second.begin(),
            it_end = it_subsection->second.end(); it != it_end; ++it)
        {
            if(whiteset.find(it->first) == whiteset.end())
            {
                std::cerr << "[IniFile] Not-whitelisted parameter \"" << it->first << "\" in section \"" << it_section->first
                          << (it_subsection->first.empty() ? "" : std::string(".") + it_subsection->first) << "\"" << std::endl;
                return false;
            }
        }
    }
    return true;
}
