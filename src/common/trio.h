#if !defined(TRIO_H)
#define TRIO_H

template<typename T1, typename T2, typename T3>
struct trio
{
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;

    T1 first;
    T2 second;
    T3 third;

    trio(const T1 & a = T1(), const T2 & b = T2(), const T3 & c = T3())
    : first(a), second(b), third(c) { }

    trio(const trio & t)
    : first(t.first), second(t.second), third(t.third) { }

    template<typename U1, typename U2, typename U3>
    trio(const trio<U1, U2, U3> & t)
    : first(t.first), second(t.second), third(t.third) { }

    trio & operator = (const trio & t)
    {
        first = t.first;
        second = t.second;
        third = t.third;
        return * this;
    }
};

template<typename T1, typename T2, typename T3>
inline bool operator == (const trio<T1, T2, T3> & x, const trio<T1, T2, T3> & y)
{
    return x.first == y.first && x.second == y.second && x.third = y.third;
}

template<typename T1, typename T2, typename T3>
inline bool operator < (const trio<T1, T2, T3> & x, const trio<T1, T2, T3> & y)
{
    return x.first < y.first ||
           (!(y.first < x.first) && x.second < y.second) ||
           (!(y.second < x.second) && x.third < y.third);
}

template<typename T1, typename T2, typename T3>
inline bool operator != (const trio<T1, T2, T3> & x, const trio<T1, T2, T3> & y)
{
    return !(x == y);
}

template<typename T1, typename T2, typename T3>
inline bool operator > (const trio<T1, T2, T3> & x, const trio<T1, T2, T3> & y)
{
    return y < x;
}

template<typename T1, typename T2, typename T3>
inline bool operator <= (const trio<T1, T2, T3> & x, const trio<T1, T2, T3> & y)
{
    return !(y < x);
}

template<typename T1, typename T2, typename T3>
inline bool operator >= (const trio<T1, T2, T3> & x, const trio<T1, T2, T3> & y)
{
    return !(x < y);
}

template<typename T1, typename T2, typename T3>
inline trio<T1, T2, T3> make_trio(const T1 & a, const T2 & b, const T3 & c)
{
    return trio<T1, T2, T3>(a, b, c);
}

#endif // TRIO_H
