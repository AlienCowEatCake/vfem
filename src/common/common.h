#if !defined COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#if defined _MSC_VER
#if !defined _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#else
#include <cmath>
#endif

#if !defined M_PI
#define M_PI 3.14159265358979323846
#endif

#if __cplusplus >= 201103L || defined __GXX_EXPERIMENTAL_CXX0X__ || (defined _MSC_VER && _MSC_VER >= 1700)
#define USE_CXX11
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cfloat>
#include <algorithm>
#include <cstring>
#include <vector>
#include <set>
#include <string>
#include <complex>
#include <map>
#include <utility>
#include <cassert>

using namespace std;

static const int IO_FILE_ERROR      = 0xfb;

#define MAYBE_UNUSED(a) (void)(a)

inline void show_progress(const char * message, size_t index, size_t num_all)
{
    MAYBE_UNUSED(message);
    MAYBE_UNUSED(index);
    MAYBE_UNUSED(num_all);
/**/
    const size_t tick = 1000;
    const char progress[4] = {'|', '/', '-', '\\'};
    index++;
    if(!(index % tick))
    {
        cout << "  ";
        if(strlen(message))
            cout << message << ": ";
        cout << index << progress[(index / tick) % 4] << num_all << '\r' << flush;
    }
    if(index >= num_all)
    {
        cout << "  ";
        if(strlen(message))
            cout << message << ": ";
        cout << num_all << progress[1] << num_all << endl;
    }
/**/
}

inline bool is_nan(double x)
{
    return x != x;
}

inline bool is_inf(double x)
{
    return !is_nan(x) && is_nan(x - x);
}

inline bool is_fpu_error(double x)
{
    return is_nan(x) || is_nan(x - x);
}

#endif // COMMON_H_INCLUDED
