#include "inifile.h"
#include <fstream>
#include <algorithm>
#include <cctype>

std::string trim(const std::string & str)
{
    size_t start = str.find_first_not_of(" \t\f\v\n\r");
    size_t stop = str.find_last_not_of(" \t\f\v\n\r");
    if(start == std::string::npos || stop == std::string::npos)
        return "";
    return str.substr(start, stop - start + 1);
}

std::string to_lowercase(const std::string & str)
{
    std::string result = str;
    std::transform(str.begin(), str.end(), result.begin(), ::tolower);
    return result;
}

inifile::inifile()
{
    status = false;
}

inifile::inifile(const std::string & filename)
{
    load(filename);
}

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
    return true;
}
