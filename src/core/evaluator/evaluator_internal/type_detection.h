#if !defined(EVALUATOR_TYPE_DETECTION_H)
#define EVALUATOR_TYPE_DETECTION_H

//#define EVALUATOR_NO_RTTI

#include <complex>
#include <string>

#if !defined(EVALUATOR_NO_RTTI)
#include <typeinfo>
#endif

namespace evaluator_internal
{

// =============================================================================================

inline bool is_floating(const float *)
{
    return true;
}

inline bool is_floating(const double *)
{
    return true;
}

inline bool is_floating(const long double *)
{
    return true;
}

template<typename T>
inline bool is_floating(const T *)
{
    return false;
}

template<typename T>
inline bool is_complex(const std::complex<T> *)
{
    return true;
}

template<typename T>
inline bool is_complex(const T *)
{
    return false;
}

template<typename T>
inline bool is_floating_complex(const std::complex<T> *)
{
    return is_floating(static_cast<T*>(NULL));
}

template<typename T>
inline bool is_floating_complex(const T *)
{
    return false;
}

// =============================================================================================

template<typename T>
inline bool is_float(const T *)
{
    return is_floating(static_cast<T*>(NULL)) && sizeof(T) == 4;
}

template<typename T>
inline bool is_double(const T *)
{
    return is_floating(static_cast<T*>(NULL)) && sizeof(T) == 8;
}

template<typename T>
inline bool is_complex_float(const std::complex<T> *)
{
    return sizeof(T) == 4;
}

template<typename T>
inline bool is_complex_float(const T *)
{
    return false;
}

template<typename T>
inline bool is_complex_double(const std::complex<T> *)
{
    return sizeof(T) == 8;
}

template<typename T>
inline bool is_complex_double(const T *)
{
    return false;
}

// =============================================================================================

#if defined(EVALUATOR_NO_RTTI)

#define ADD_GET_TYPE_NAME(TYPE) \
    inline std::string get_type_name(const TYPE *) \
    { \
        return std::string(#TYPE); \
    }
    ADD_GET_TYPE_NAME(signed char)
    ADD_GET_TYPE_NAME(unsigned char)
    ADD_GET_TYPE_NAME(char)
    ADD_GET_TYPE_NAME(signed short)
    ADD_GET_TYPE_NAME(unsigned short)
    ADD_GET_TYPE_NAME(signed int)
    ADD_GET_TYPE_NAME(unsigned int)
    ADD_GET_TYPE_NAME(signed long)
    ADD_GET_TYPE_NAME(unsigned long)
    ADD_GET_TYPE_NAME(signed long long) /// @note C++11
    ADD_GET_TYPE_NAME(unsigned long long) /// @note C++11
    ADD_GET_TYPE_NAME(float)
    ADD_GET_TYPE_NAME(double)
    ADD_GET_TYPE_NAME(long double)
    ADD_GET_TYPE_NAME(std::complex<float>)
    ADD_GET_TYPE_NAME(std::complex<double>)
    ADD_GET_TYPE_NAME(std::complex<long double>)
#undef ADD_GET_TYPE_NAME

template<typename T>
inline std::string get_type_name(const T *)
{
    return "unknown";
}

#else

template<typename T>
inline std::string get_type_name(const T *)
{
    return typeid(T).name();
}

#endif

// =============================================================================================

} // namespace evaluator_internal

#endif // EVALUATOR_TYPE_DETECTION_H

