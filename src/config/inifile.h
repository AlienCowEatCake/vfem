#ifndef INIFILE_H
#define INIFILE_H

#include <map>
#include <string>
#include <sstream>
#include <list>

std::string trim(const std::string & str);
std::string to_lowercase(const std::string & str);

class inifile
{
public:
    inifile();
    inifile(const std::string & filename);
    bool load(const std::string & filename);
    inline bool good() const
    {
        return status;
    }

    template<typename U>
    bool get(const std::string & section, const U & subsection, const std::string & parameter, const bool & fallback) const
    {
        std::string subsection_str;
        std::stringstream subsection_sst;
        subsection_sst << subsection;
        subsection_sst >> subsection_str;
        std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::const_iterator it_section = values.find(section);
        if(it_section != values.end())
        {
            std::map<std::string, std::map<std::string, std::string> >::const_iterator it_subsection = it_section->second.find(subsection_str);
            if(it_subsection != it_section->second.end())
            {
                std::map<std::string, std::string>::const_iterator it_parameter = it_subsection->second.find(parameter);
                if(it_parameter != it_subsection->second.end())
                {
                    std::string value = to_lowercase(it_parameter->second);
                    if(value == "yes" || value == "true" || value == "1")
                        return true;
                    else
                        return false;
                }
            }
        }
        return fallback;
    }

    template<typename T, typename U>
    T get(const std::string & section, const U & subsection, const std::string & parameter, const T & fallback) const
    {
        std::string subsection_str;
        std::stringstream subsection_sst;
        subsection_sst << subsection;
        subsection_sst >> subsection_str;
        std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::const_iterator it_section = values.find(section);
        if(it_section != values.end())
        {
            std::map<std::string, std::map<std::string, std::string> >::const_iterator it_subsection = it_section->second.find(subsection_str);
            if(it_subsection != it_section->second.end())
            {
                std::map<std::string, std::string>::const_iterator it_parameter = it_subsection->second.find(parameter);
                if(it_parameter != it_subsection->second.end())
                {
                    std::stringstream sst(it_parameter->second);
                    T result;
                    sst >> result;
                    return result;
                }
            }
        }
        return fallback;
    }

    template<typename T>
    std::list<T> enumerate(const std::string & section, const T & type) const
    {
        (void)(type);
        std::list<std::string> result;
        typename std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::const_iterator it_section = values.find(section);
        if(it_section != values.end())
        {
            for(std::map<std::string, std::map<std::string, std::string> >::const_iterator it = it_section->second.begin(); it != it_section->second.end(); ++it)
            {
                std::stringstream sst(it->first);
                T subsection = T();
                sst >> subsection;
                result.push_back(subsection);
            }
        }
        return result;
    }

protected:

    std::map<std::string, std::map<std::string, std::map<std::string, std::string> > > values;

    bool status;
};


#endif // INIFILE_H
