#if !defined(UTILS_CXXVERSION_H_INCLUDED)
#define UTILS_CXXVERSION_H_INCLUDED

/**
  * @def USE_CXX11 будет определен только для тех компиляторов, которые умеют C++11
  */
#if __cplusplus >= 201103L || \
    (defined(__GNUC__) && defined(__GNUC_MINOR__) && defined(__GXX_EXPERIMENTAL_CXX0X__) && \
    (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8))) || \
    (defined(_MSC_VER) && _MSC_VER >= 1800)
#if !defined(__APPLE__) || (defined(__APPLE__) && defined(MAC_OS_X_VERSION_MAX_ALLOWED) && \
    (MAC_OS_X_VERSION_MAX_ALLOWED > MAC_OS_X_VERSION_10_6))
#define USE_CXX11
#endif
#endif

#endif // UTILS_CXXVERSION_H_INCLUDED
