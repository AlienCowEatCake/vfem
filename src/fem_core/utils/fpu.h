#if !defined(UTILS_FPU_H_INCLUDED)
#define UTILS_FPU_H_INCLUDED

#include <cmath>

#include "cxxversion.h"

namespace fem_core { namespace utils { namespace fpu {

/**
 * @brief Determines if the given floating point number arg is a not-a-number (NaN) value.
 * @param[in] x floating point value.
 * @return true if arg is a NaN, false otherwise.
 */
template<typename T>
inline bool is_nan(const T & x)
{
#if defined(USE_CXX11)
    return std::isnan(x);
#else
    return x != x;
#endif
}

/**
 * @brief Determines if the given floating point number arg is a positive or negative infinity.
 * @param[in] x floating point value.
 * @return true if arg is infinite, false otherwise.
 */
template<typename T>
inline bool is_inf(const T & x)
{
#if defined(USE_CXX11)
    return std::isinf(x);
#else
    return !is_nan(x) && is_nan(x - x);
#endif
}

/**
 * @brief Determines if the given floating point number arg is a not-a-number (NaN) or infinity value.
 * @param[in] x floating point value.
 * @return true if arg is a NaN or infinite, false otherwise.
 */
template<typename T>
inline bool is_fpu_error(const T & x)
{
#if defined(USE_CXX11)
    return !std::isfinite(x);
#else
    return is_nan(x) || is_nan(x - x);
#endif
}

}}} // namespace fem_core::utils::fpu

#endif // UTILS_FPU_H_INCLUDED
