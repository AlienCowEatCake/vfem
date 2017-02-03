#if !defined(CONTAINERS_GENERIC_TRIO_H_INCLUDED)
#define CONTAINERS_GENERIC_TRIO_H_INCLUDED

namespace fem_core { namespace containers { namespace generic {

/**
 * @brief The trio struct is a struct template that provides a way to store three heterogeneous objects as a single unit.
 */
template<typename T1, typename T2, typename T3>
struct trio
{
    /**
     * @brief Type of member first.
     */
    typedef T1 first_type;
    /**
     * @brief Type of member second.
     */
    typedef T2 second_type;
    /**
     * @brief Type of member third.
     */
    typedef T3 third_type;

    /**
     * @brief The first value in the trio.
     */
    T1 first;
    /**
     * @brief The second value in the trio.
     */
    T2 second;
    /**
     * @brief The third value in the trio.
     */
    T3 third;

    /**
     * @brief Constructs a trio object.
     * @param[in] a An object of the type of first, or some other type implicitly convertible to it.
     * @param[in] b An object of the type of second, or some other type implicitly convertible to it.
     * @param[in] c An object of the type of third, or some other type implicitly convertible to it.
     */
    trio(const T1 & a = T1(), const T2 & b = T2(), const T3 & c = T3())
    : first(a), second(b), third(c) { }

    /**
     * @brief The object is initialized with the contents of the t trio object.
     * @param[in] t Another trio object.
     */
    trio(const trio & t)
    : first(t.first), second(t.second), third(t.third) { }

    /**
     * @brief The object is initialized with the contents of the t trio object.
     * @param[in] t Another trio object with some other type implicitly convertible to it.
     */
    template<typename U1, typename U2, typename U3>
    trio(const trio<U1, U2, U3> & t)
    : first(t.first), second(t.second), third(t.third) { }

    /**
     * @brief Assigns t as the new content for the trio object.
     * @param[in] t Another trio object.
     * @return *this.
     */
    trio & operator = (const trio & t)
    {
        if(this != & t)
        {
            first = t.first;
            second = t.second;
            third = t.third;
        }
        return * this;
    }

    /**
     * @brief Assigns t as the new content for the trio object.
     * @param[in] t Another trio object with some other type implicitly convertible to it.
     * @return *this.
     */
    template<typename U1, typename U2, typename U3>
    trio & operator = (const trio<U1, U2, U3> & t)
    {
        if(this != & t)
        {
            first = t.first;
            second = t.second;
            third = t.third;
        }
        return * this;
    }
};

/**
 * @brief Checks if the values of two operands are equal or not.
 * @param[in] x Left trio operand.
 * @param[in] y Right trio operand.
 * @return true or false.
 */
template<typename T1, typename T2, typename T3>
inline bool operator == (const trio<T1, T2, T3> & x, const trio<T1, T2, T3> & y)
{
    return x.first == y.first && x.second == y.second && x.third = y.third;
}

/**
 * @brief Checks if the value of left operand is less than the value of right operand.
 * @param[in] x Left trio operand.
 * @param[in] y Right trio operand.
 * @return true or false.
 */
template<typename T1, typename T2, typename T3>
inline bool operator < (const trio<T1, T2, T3> & x, const trio<T1, T2, T3> & y)
{
    return x.first < y.first ||
           (!(y.first < x.first) && x.second < y.second) ||
           (!(y.second < x.second) && x.third < y.third);
}

/**
 * @brief Checks if the values of two operands are equal or not.
 * @param[in] x Left trio operand.
 * @param[in] y Right trio operand.
 * @return true or false.
 */
template<typename T1, typename T2, typename T3>
inline bool operator != (const trio<T1, T2, T3> & x, const trio<T1, T2, T3> & y)
{
    return !(x == y);
}

/**
 * @brief Checks if the value of left operand is greater than the value of right operand.
 * @param[in] x Left trio operand.
 * @param[in] y Right trio operand.
 * @return true or false.
 */
template<typename T1, typename T2, typename T3>
inline bool operator > (const trio<T1, T2, T3> & x, const trio<T1, T2, T3> & y)
{
    return y < x;
}

/**
 * @brief Checks if the value of left operand is less than or equal to the value of right operand.
 * @param[in] x Left trio operand.
 * @param[in] y Right trio operand.
 * @return true or false.
 */
template<typename T1, typename T2, typename T3>
inline bool operator <= (const trio<T1, T2, T3> & x, const trio<T1, T2, T3> & y)
{
    return !(y < x);
}

/**
 * @brief Checks if the value of left operand is greater than or equal to the value of right operand.
 * @param[in] x Left trio operand.
 * @param[in] y Right trio operand.
 * @return true or false.
 */
template<typename T1, typename T2, typename T3>
inline bool operator >= (const trio<T1, T2, T3> & x, const trio<T1, T2, T3> & y)
{
    return !(x < y);
}

/**
 * @brief Constructs a trio object with its first element set to a, its second element set to b, and its third element set to c.
 * @param[in] a Value for the members first of the trio object being constructed.
 * @param[in] b Value for the second first of the trio object being constructed.
 * @param[in] c Value for the third first of the trio object being constructed.
 * @return A trio object.
 */
template<typename T1, typename T2, typename T3>
inline trio<T1, T2, T3> make_trio(const T1 & a, const T2 & b, const T3 & c)
{
    return trio<T1, T2, T3>(a, b, c);
}

}}} // namespace fem_core::containers::generic

#endif // CONTAINERS_GENERIC_TRIO_H_INCLUDED
