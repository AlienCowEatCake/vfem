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

#if defined _MSC_VER && _MSC_VER <= 1200
#define for  if (0) {} else for
#endif

#if !defined M_PI
#define M_PI 3.14159265358979323846
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

static const int ADDRESSING_ERROR   = 0xfa;
static const int IO_FILE_ERROR      = 0xfb;
static const int NULL_PTR_ERROR     = 0xfc;

#define MAYBE_UNUSED(a) (void)(a)

inline void show_progress(const char * message, size_t index, size_t num_all)
{
    MAYBE_UNUSED(message);
    MAYBE_UNUSED(index);
    MAYBE_UNUSED(num_all);
/**/
    const size_t tick = 100;
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

#if defined __WATCOMC__

istream & getline(istream & is, string & str)
{
    const int bufsize = 255;
    char buf[bufsize];
    is.getline(buf, bufsize);
    str.assign(buf);
    return is;
}

ostream & operator << (ostream & os, const string & str)
{
    os << str.c_str();
    return os;
}

ostream & operator << (ostream & os, const complex<double> & c)
{
    os << '(' << c.real() << ',' << c.imag() << ')';
    return os;
}
#endif

#endif // COMMON_H_INCLUDED
