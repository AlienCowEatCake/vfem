#include "inifile.h"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <set>

namespace core { namespace utils {

/**
 * @brief Загрузка ini-файла
 * @param[in] filename имя входного файла
 * @return true - файл успешно считан, false - возникли ошибки
 */
bool inifile::load(const std::string & filename)
{
    std::ifstream ifs(filename.c_str());
    if(!ifs.good())
    {
        m_status = false;
        return false;
    }

    std::string line;
    std::getline(ifs, line);
    line = strings::trim(line);
    do
    {
        if(line.length() > 1 && line[0] == '[')
        {
            line = strings::to_lowercase(line.substr(1, line.length() - 2));
            std::string section = line, subsection;
            size_t dot_pos = line.find_first_of(".");
            if(dot_pos != std::string::npos)
            {
                section = line.substr(0, dot_pos);
                subsection = line.substr(dot_pos + 1);
            }
            m_values[section]._[subsection]; ///< Пустая секция - тоже валидная секция

            do
            {
                getline(ifs, line);
                line = strings::trim(line);
                if(line.length() > 1 && line[0] != ';' && line[0] != '#')
                {
                    size_t eq_pos = line.find_first_of("=");
                    if(eq_pos != std::string::npos)
                    {
                        std::string param = strings::to_lowercase(strings::trim(line.substr(0, eq_pos)));
                        std::string value = strings::trim(line.substr(eq_pos + 1));
                        if(value.length() > 1 && value[0] == '\"')
                            value = strings::trim(value.substr(1, value.length() - 2));
                        m_values[section]._[subsection]._[param] = value;
                    }
                }
            }
            while(ifs.good() && !(line.length() > 1 && line[0] == '['));
        }
        else
        {
            getline(ifs, line);
            line = strings::trim(line);
        }
    }
    while(ifs.good());

    ifs.close();
    m_status = true;
    return true;
}

/**
 * @brief Получить значение параметра с именем "parameter" в секции "section" подсекции "subsection"
 * @param[in] section секция
 * @param[in] subsection подсекция
 * @param[in] parameter параметр
 * @return значение параметра или NULL, если такого параметра нет
 */
const std::string * inifile::get_internal(const std::string & section, const std::string & subsection, const std::string & parameter) const
{
    // Найдем искомую секцию
    std::map<std::string, section_t>::const_iterator it_section = m_values.find(strings::to_lowercase(section));
    if(it_section != m_values.end())
    {
        // Найдем искомую подсекцию
        std::map<std::string, subsection_t>::const_iterator it_subsection = it_section->second._.find(strings::to_lowercase(subsection));
        if(it_subsection != it_section->second._.end())
        {
            // И найдем искомый параметр
            std::map<std::string, std::string>::const_iterator it_parameter = it_subsection->second._.find(strings::to_lowercase(parameter));
            if(it_parameter != it_subsection->second._.end())
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
 * @param[in] in входной std::list<std::string>
 * @return выходной std::set<std::string>
 */
std::set<std::string> list_to_set_and_lower(const std::list<std::string> & in)
{
    std::set<std::string> out;
    for(std::list<std::string>::const_iterator it = in.begin(), it_end = in.end(); it != it_end; ++it)
        out.insert(strings::to_lowercase(* it));
    return out;
}

}

/**
 * @brief Проверить, что в файле нет никаких секций, кроме секций из списка "whitelist"
 * @param[in] whitelist список разрешенных секций
 * @return true - все секции присутствуют в "whitelist", false - есть неразрешенная секция
 * @note Полезно для контроля опечаток, так как в остальных местах это штатная ситуация
 */
bool inifile::check_sections(const std::list<std::string> & whitelist) const
{
    if(!m_status)
        return false;
    std::set<std::string> whiteset = list_to_set_and_lower(whitelist);
    for(std::map<std::string, section_t>::const_iterator
        it = m_values.begin(), it_end = m_values.end(); it != it_end; ++it)
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
 * @param[in] section секция, для которой будет выполнена проверка
 * @param[in] whitelist список разрешенных параметров
 * @return true - все параметры присутствуют в "whitelist", false - есть неразрешенный параметр
 * @note Полезно для контроля опечаток, так как в остальных местах это штатная ситуация
 */
bool inifile::check_parameters(const std::string & section, const std::list<std::string> & whitelist) const
{
    if(!m_status)
        return false;
    std::map<std::string, section_t>::const_iterator it_section = m_values.find(strings::to_lowercase(section));
    if(it_section == m_values.end())
        return true;
    std::set<std::string> whiteset = list_to_set_and_lower(whitelist);
    for(std::map<std::string, subsection_t>::const_iterator it_subsection = it_section->second._.begin(),
        it_subsection_end = it_section->second._.end(); it_subsection != it_subsection_end; ++it_subsection)
    {
        for(std::map<std::string, std::string>::const_iterator it = it_subsection->second._.begin(),
            it_end = it_subsection->second._.end(); it != it_end; ++it)
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

}} // namespace core::utils
