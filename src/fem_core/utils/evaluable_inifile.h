#if !defined(UTILS_EVALUABLE_INIFILE_H_INCLUDED)
#define UTILS_EVALUABLE_INIFILE_H_INCLUDED

#include <complex>
#include "inifile.h"
#include "../evaluator/evaluator.h"

namespace fem_core { namespace utils {

/**
 * @brief Класс для работы с ini-файлами с поддержкой вычисляемых выражений
 */
class evaluable_inifile : public inifile
{
public:
    evaluable_inifile() : inifile() {}

    evaluable_inifile(const std::string & filename) : inifile(filename) {}

#define ADD_GET_FOR_TYPE(TYPE) \
    template<typename U> \
    TYPE get(const std::string & section, const U & subsection, const std::string & parameter, const TYPE & fallback) const \
    { \
        return evaluate_and_get(section, subsection, parameter, fallback); \
    }
    ADD_GET_FOR_TYPE(float)
    ADD_GET_FOR_TYPE(double)
    ADD_GET_FOR_TYPE(long double)
    ADD_GET_FOR_TYPE(signed short)
    ADD_GET_FOR_TYPE(unsigned short)
    ADD_GET_FOR_TYPE(signed int)
    ADD_GET_FOR_TYPE(unsigned int)
    ADD_GET_FOR_TYPE(signed long)
    ADD_GET_FOR_TYPE(unsigned long)
    ADD_GET_FOR_TYPE(signed long long)
    ADD_GET_FOR_TYPE(unsigned long long)
    ADD_GET_FOR_TYPE(std::complex<float>)
    ADD_GET_FOR_TYPE(std::complex<double>)
    ADD_GET_FOR_TYPE(std::complex<long double>)
#undef ADD_GET_FOR_TYPE

    template<typename U>
    inline std::string get(const std::string & section, const U & subsection, const std::string & parameter, const char * fallback) const
    {
        return inifile::get(section, subsection, parameter, fallback);
    }

    template<typename T, typename U>
    T get(const std::string & section, const U & subsection, const std::string & parameter, const T & fallback) const
    {
        return inifile::get(section, subsection, parameter, fallback);
    }

protected:

    template<typename T, typename U>
    T evaluate_and_get(const std::string & section, const U & subsection, const std::string & parameter, const T & fallback) const
    {
        const std::string * str = get_internal(section, subsection_to_str(subsection), parameter);
        if(!str)
            return fallback;
        evaluator<T> eval(*str);
        T result;
        if(eval.is_parsed() && eval.calculate(result))
            return result;
        return fallback;
    }

};

}} // namespace fem_core::utils

#endif // UTILS_EVALUABLE_INIFILE_H_INCLUDED
