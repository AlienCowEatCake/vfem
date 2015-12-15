#ifndef PARSER_H
#define PARSER_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <stack>
#include <deque>
#include <algorithm>
#include <cmath>
#include <utility>
#include <complex>
#include <typeinfo>
#include <limits>
#include <cctype>

#if defined(_MSC_VER) && _MSC_VER < 1800
namespace std
{
    // http://functions.wolfram.com/ElementaryFunctions/ArcCosh/02/
    template<typename T>
    T acosh(T x) { return log(x + sqrt(x * x - 1)); }
    // http://functions.wolfram.com/ElementaryFunctions/ArcSinh/02/
    template<typename T>
    T asinh(T x) { return log(x + sqrt(x * x + 1)); }
    // http://functions.wolfram.com/ElementaryFunctions/ArcTanh/02/
    template<typename T>
    T atanh(T x) { return  0.5 * log((1.0 + x) / (1.0 - x)); }
}
#endif

namespace parser_internal
{
    using namespace std;

    template<typename T> T pi_sin   (const T & arg) { return sin(arg);   }
    template<typename T> T pi_cos   (const T & arg) { return cos(arg);   }
    template<typename T> T pi_tan   (const T & arg) { return tan(arg);   }
    template<typename T> T pi_sinh  (const T & arg) { return sinh(arg);  }
    template<typename T> T pi_cosh  (const T & arg) { return cosh(arg);  }
    template<typename T> T pi_tanh  (const T & arg) { return tanh(arg);  }
    template<typename T> T pi_log   (const T & arg) { return log(arg);   }
    template<typename T> T pi_log10 (const T & arg) { return log10(arg); }
    template<typename T> T pi_exp   (const T & arg) { return exp(arg);   }
    template<typename T> T pi_sqrt  (const T & arg) { return sqrt(arg);  }
    template<typename T> complex<T> pi_abs   (const complex<T> & arg) { return abs(arg);          }
    template<typename T> complex<T> pi_asin  (const complex<T> & arg) { return asin(arg.real());  }
    template<typename T> complex<T> pi_acos  (const complex<T> & arg) { return acos(arg.real());  }
    template<typename T> complex<T> pi_atan  (const complex<T> & arg) { return atan(arg.real());  }
    template<typename T> complex<T> pi_asinh (const complex<T> & arg) { return asinh(arg.real()); }
    template<typename T> complex<T> pi_acosh (const complex<T> & arg) { return acosh(arg.real()); }
    template<typename T> complex<T> pi_atanh (const complex<T> & arg) { return atanh(arg.real()); }
    template<typename T> complex<T> pi_imag  (const complex<T> & arg) { return arg.imag();        }
    template<typename T> complex<T> pi_real  (const complex<T> & arg) { return arg.real();        }
    template<typename T> complex<T> pi_conj  (const complex<T> & arg) { return conj(arg);         }
    template<typename T> T pi_abs   (const T & arg) { return fabs(arg);  }
    template<typename T> T pi_asin  (const T & arg) { return asin(arg);  }
    template<typename T> T pi_acos  (const T & arg) { return acos(arg);  }
    template<typename T> T pi_atan  (const T & arg) { return atan(arg);  }
    template<typename T> T pi_asinh (const T & arg) { return asinh(arg); }
    template<typename T> T pi_acosh (const T & arg) { return acosh(arg); }
    template<typename T> T pi_atanh (const T & arg) { return atanh(arg); }
    template<typename T> T pi_imag  (const T & arg) { return arg * (T)0; }
    template<typename T> T pi_real  (const T & arg) { return arg;        }
    template<typename T> T pi_conj  (const T & arg) { return arg;        }

    template<typename T>
    void init_functions(map<string, T(*)(const T &)> & funcs_map)
    {
        funcs_map["imag"]  = pi_imag;
        funcs_map["real"]  = pi_real;
        funcs_map["conj"]  = pi_conj;
        funcs_map["sin"]   = pi_sin;
        funcs_map["cos"]   = pi_cos;
        funcs_map["tan"]   = pi_tan;
        funcs_map["asin"]  = pi_asin;
        funcs_map["acos"]  = pi_acos;
        funcs_map["atan"]  = pi_atan;
        funcs_map["sinh"]  = pi_sinh;
        funcs_map["cosh"]  = pi_cosh;
        funcs_map["tanh"]  = pi_tanh;
        funcs_map["asinh"] = pi_asinh;
        funcs_map["acosh"] = pi_acosh;
        funcs_map["atanh"] = pi_atanh;
        funcs_map["log"]   = pi_log;
        funcs_map["log10"] = pi_log10;
        funcs_map["abs"]   = pi_abs;
        funcs_map["exp"]   = pi_exp;
        funcs_map["sqrt"]  = pi_sqrt;
    }

    template<typename T> T pi_plus  (const T & larg, const T & rarg) { return larg + rarg; }
    template<typename T> T pi_minus (const T & larg, const T & rarg) { return larg - rarg; }
    template<typename T> T pi_mult  (const T & larg, const T & rarg) { return larg * rarg; }
    template<typename T> T pi_div   (const T & larg, const T & rarg) { return larg / rarg; }
    template<typename T> T pi_pow   (const T & larg, const T & rarg) { return pow(larg, rarg); }

    template<typename T>
    void init_operators(map<char, pair<unsigned short int, T(*)(const T &, const T &)> > & opers_map)
    {
        typedef pair<unsigned short int, T(*)(const T &, const T &)> oper_type;
        opers_map['+'] = oper_type(1, pi_plus);
        opers_map['-'] = oper_type(1, pi_minus);
        opers_map['*'] = oper_type(2, pi_mult);
        opers_map['/'] = oper_type(2, pi_div);
        opers_map['^'] = oper_type(3, pi_pow);
    }

    template<typename T>
    void init_constants(map<string, T> & consts_map)
    {
        consts_map["pi"] = static_cast<T>(3.14159265358979323846264338327950);
        consts_map["e"]  = static_cast<T>(2.71828182845904523536028747135266);
        if(typeid(T) == typeid(complex<float>) ||
           typeid(T) == typeid(complex<double>) ||
           typeid(T) == typeid(complex<long double>))
        {
            consts_map["i"] = sqrt(static_cast<T>(-1.0));
            consts_map["j"] = sqrt(static_cast<T>(-1.0));
        }
    }

    template<typename T>
    class parser_object
    {
    protected:
        enum parser_object_type
        {
            PI_OBJ_OPERATOR,
            PI_OBJ_FUNCTION,
            PI_OBJ_VARIABLE,
            PI_OBJ_CONSTANT
        };
        parser_object_type type;
        string str_;
        T(* func)(const T &);
        T(* oper)(const T &, const T &);
        T value;
        void init(parser_object_type n_type, const string & n_str,
                  T(* n_func)(const T &), T(* n_oper)(const T &, const T &), const T & n_value)
        {
            type = n_type;
            str_ = n_str;
            value = n_value;
            func = n_func;
            oper = n_oper;
        }
    public:
        inline bool is_variable() const { return type == PI_OBJ_VARIABLE; }
        inline bool is_constant() const { return type == PI_OBJ_CONSTANT; }
        inline bool is_function() const { return type == PI_OBJ_FUNCTION; }
        inline bool is_operator() const { return type == PI_OBJ_OPERATOR; }
        inline const string & str() const { return str_; }
        inline T eval() const { return value; }
        inline T eval(const T & arg) const { return func(arg); }
        inline T eval(const T & larg, const T & rarg) const { return oper(larg, rarg); }
        parser_object(const string & new_str)
        {
            init(PI_OBJ_VARIABLE, new_str, NULL, NULL, T());
        }
        parser_object(const string & new_str, const T & new_value)
        {
            init(PI_OBJ_CONSTANT, new_str, NULL, NULL, new_value);
        }
        parser_object(const string & new_str, T(* new_func)(const T &))
        {
            init(PI_OBJ_FUNCTION, new_str, new_func, NULL, T());
        }
        parser_object(const string & new_str, T(* new_oper)(const T &, const T &))
        {
            init(PI_OBJ_OPERATOR, new_str, NULL, new_oper, T());
        }
    };
}

template<typename T>
class parser
{
protected:
    std::vector<std::string> expression;
    std::vector<parser_internal::parser_object<T> > expression_objects;
    std::map<std::string, T(*)(const T &)> functions;
    std::map<std::string, T> constants;
    std::map<char, std::pair<unsigned short int, T(*)(const T &, const T &)> > operators;
    bool status;
    std::string error_string;

    void init()
    {
        using namespace std;
        using namespace parser_internal;
        init_functions(functions);
        init_operators(operators);
        init_constants(constants);
    }

    void convert_to_objects()
    {
        using namespace std;
        using namespace parser_internal;
        expression_objects.clear();
        expression_objects.reserve(expression.size());
        for(vector<string>::const_iterator it = expression.begin(); it != expression.end(); ++it)
        {
            typename map<string, T>::const_iterator itc = constants.find(*it);
            if(itc != constants.end())
            {
                expression_objects.push_back(parser_object<T>(*it, itc->second));
            }
            else
            {
                typename map<char, pair<unsigned short int, T(*)(const T &, const T &)> >::const_iterator ito = operators.find((*it)[0]);
                if(ito != operators.end())
                {
                    expression_objects.push_back(parser_object<T>(*it, ito->second.second));
                }
                else
                {
                    typename map<string, T(*)(const T &)>::const_iterator itf = functions.find(*it);
                    if(itf != functions.end())
                    {
                        expression_objects.push_back(parser_object<T>(*it, itf->second));
                    }
                    else
                    {
                        expression_objects.push_back(parser_object<T>(*it));
                    }
                }
            }
        }
        expression.clear();
    }

public:
    parser()
    {
        init();
        status = false;
    }

    parser(const std::string & str)
    {
        init();
        parse(str);
    }

    const std::string & get_error() const
    {
        return error_string;
    }

    void set_const(const std::string & name, const T & value)
    {
        constants[name] = value;
    }

    void reset_const()
    {
        using namespace parser_internal;
        constants.clear();
        init_constants(constants);
    }

    bool is_parsed() const
    {
        return status;
    }

    bool parse(const std::string & str)
    {
        using namespace std;
        using namespace parser_internal;

        expression.clear();
        expression_objects.clear();
        error_string = "";

        status = true;
        bool str_begin = true;
        bool unary_minus = false;
        stack<string> st;

        for(string::const_iterator it = str.begin(); it != str.end() && status;)
        {
            char sym = *it;
            if(sym >= '0' && sym <= '9')
            {
                string a;
                while(it != str.end() && sym >= '0' && sym <= '9')
                {
                    a.push_back(sym);
                    ++it;
                    if(it != str.end())
                    {
                        sym = *it;
                    }
                }
                if(it != str.end() && (sym == '.' || sym == ','))
                {
                    a.push_back('.');
                    ++it;
                    if(it != str.end())
                    {
                        sym = *it;
                    }
                    while(it != str.end() && sym >= '0' && sym <= '9')
                    {
                        a.push_back(sym);
                        ++it;
                        if(it != str.end())
                        {
                            sym = *it;
                        }
                    }
                }
                if(it != str.end() && (sym == 'e' || sym == 'E' || sym == 'd' || sym == 'D'))
                {
                    a.push_back('e');
                    ++it;
                    if(it != str.end())
                    {
                        sym = *it;
                    }
                    if(it != str.end() && (sym == '-' || sym == '+'))
                    {
                        a.push_back(sym);
                        ++it;
                        if(it != str.end())
                        {
                            sym = *it;
                        }
                    }
                    while(it != str.end() && sym >= '0' && sym <= '9')
                    {
                        a.push_back(sym);
                        ++it;
                        if(it != str.end())
                        {
                            sym = *it;
                        }
                    }
                }

                if(unary_minus)
                {
                    a = "-" + a;
                    unary_minus = false;
                }
                stringstream b;
                b << a;
                double c;
                b >> c;
                expression.push_back(a);
                constants[a] = static_cast<T>(c);
                str_begin = false;
            }
            else if(sym == '(')
            {
                st.push("(");
                str_begin = true;
                ++it;
            }
            else if(sym == ')')
            {
                while(!st.empty() && st.top() != "(")
                {
                    expression.push_back(st.top());
                    st.pop();
                }
                if(st.empty())
                {
                    status = false;
                    error_string = "Wrong brackets balance!";
                    break;
                }
                if(str_begin)
                {
                    status = false;
                    error_string = "Unexpected ')'!";
                    break;
                }
                st.pop();
                if(!st.empty() && functions.find(st.top()) != functions.end())
                {
                    expression.push_back(st.top());
                    st.pop();
                }
                str_begin = false;
                ++it;
            }
            else if(operators.find(sym) != operators.end())
            {
                if(sym == '-' && str_begin)
                {
                    string::const_iterator it2 = it + 1;
                    if(it2 == str.end())
                    {
                        status = false;
                        error_string = "Unexpected '-'!";
                        break;
                    }
                    if(*it2 >= '0' && *it2 <= '9')
                    {
                        unary_minus = true;
                    }
                    else
                    {
                        constants["-1.0"] = static_cast<T>(-1.0);
                        expression.push_back("-1.0");
                        st.push("*");
                    }
                }
                else
                {
                    char op;
                    if(!st.empty()) op = st.top().c_str()[0];
                    while(!st.empty() && operators.find(op) != operators.end() &&
                          operators[sym].first <= operators[op].first)
                    {
                        expression.push_back(st.top());
                        st.pop();
                        if(!st.empty()) op = st.top().c_str()[0];
                    }
                    string tmp;
                    tmp.push_back(sym);
                    st.push(tmp);
                }
                str_begin = false;
                ++it;
            }
            else if(sym != ' ' && sym != '\t' && sym != '\0' && sym != '\r' && sym != '\n')
            {
                string funcname;
                while(it != str.end() && *it != '(' && *it != ')' &&
                      *it != ' ' && *it != '\t' && *it != '\0' && *it != '\r' && *it != '\n' &&
                      operators.find(*it) == operators.end())
                {
                    funcname.push_back(*it);
                    ++it;
                }
                transform(funcname.begin(), funcname.end(), funcname.begin(), ::tolower);
                if(functions.find(funcname) != functions.end())
                {
                    st.push(funcname);
                }
                else if(*it != '(')
                {
                    expression.push_back(funcname);
                }
                else
                {
                    status = false;
                    error_string = "Wrong function!";
                    break;
                }
                str_begin = false;
            }
            else
            {
                ++it;
            }
        }

        while(status && !st.empty())
        {
            if(operators.find(st.top().c_str()[0]) == operators.end())
            {
                status = false;
                error_string = "Wrong expression!";
                break;
            }
            expression.push_back(st.top());
            st.pop();
        }

        if(status) convert_to_objects();
        return status;
    }

    bool simplify()
    {
        using namespace std;
        using namespace parser_internal;

        if(!is_parsed())
        {
            error_string = "Not parsed!";
            return false;
        }

        bool was_changed;
        do
        {
            deque<parser_object<T> > dq;
            was_changed = false;

            for(typename vector<parser_object<T> >::iterator it = expression_objects.begin(); it != expression_objects.end(); ++it)
            {
                if(it->is_operator())
                {
                    parser_object<T> arg2 = dq.back();
                    if(arg2.is_constant() || (arg2.is_variable() && (constants.find(arg2.str()) != constants.end())))
                    {
                        dq.pop_back();
                        parser_object<T> arg1 = dq.back();
                        if(arg1.is_constant() || (arg1.is_variable() && (constants.find(arg1.str()) != constants.end())))
                        {
                            dq.pop_back();
                            T varg1 = (arg1.is_constant() ? arg1.eval() : constants[arg1.str()]);
                            T varg2 = (arg2.is_constant() ? arg2.eval() : constants[arg2.str()]);
                            T val = it->eval(varg1, varg2);
                            stringstream sst;
                            sst.precision(17);
                            sst.setf(ios::scientific);
                            sst << val;
                            string sst_st = sst.str();
                            constants[sst_st] = val;
                            dq.push_back(parser_object<T>(sst_st, val));
                        }
                        else
                        {
                            dq.push_back(arg2);
                            dq.push_back(*it);
                        }
                    }
                    else
                    {
                        dq.push_back(*it);
                    }
                }
                else if(it->is_function())
                {
                    parser_object<T> arg = dq.back();
                    if(arg.is_constant() || (arg.is_variable() && (constants.find(arg.str()) != constants.end())))
                    {
                        dq.pop_back();
                        T varg = (arg.is_constant() ? arg.eval() : constants[arg.str()]);
                        T val = it->eval(varg);
                        stringstream sst;
                        sst.precision(17);
                        sst.setf(ios::scientific);
                        sst << val;
                        string sst_st = sst.str();
                        constants[sst_st] = val;
                        dq.push_back(parser_object<T>(sst_st, val));
                    }
                    else
                    {
                        dq.push_back(*it);
                    }
                }
                else
                {
                    dq.push_back(*it);
                }
            }

            if(expression_objects.size() > dq.size())
            {
                expression_objects.clear();
                expression_objects.reserve(dq.size());
                while(!dq.empty())
                {
                    expression_objects.push_back(dq.front());
                    dq.pop_front();
                }
                was_changed = true;
            }
            else
            {
                dq.clear();
            }

            for(typename vector<parser_object<T> >::iterator it = expression_objects.begin(); it != expression_objects.end(); ++it)
            {
                if(it->is_operator())
                {
                    parser_object<T> arg2 = dq.back();
                    dq.pop_back();
                    parser_object<T> arg1 = dq.back();
                    // Such things as a*0 or 0*a
                    if(it->str() == "*" && ((arg2.is_constant() && arg2.eval() == static_cast<T>(0.0)) ||
                                            (arg1.is_constant() && arg1.eval() == static_cast<T>(0.0) &&
                                             !arg2.is_operator() && !arg2.is_function())))
                    {
                        dq.pop_back();
                        if(arg2.is_constant() && arg2.eval() == static_cast<T>(0.0))
                            dq.push_back(arg2);
                        else
                            dq.push_back(arg1);
                    }
                    // Such things as a*1 or 1*a
                    else if(it->str() == "*" && ((arg2.is_constant() && arg2.eval() == static_cast<T>(1.0)) ||
                                                 (arg1.is_constant() && arg1.eval() == static_cast<T>(1.0) &&
                                                  !arg2.is_operator() && !arg2.is_function())))
                    {
                        dq.pop_back();
                        if(arg2.is_constant() && arg2.eval() == static_cast<T>(1.0))
                            dq.push_back(arg1);
                        else
                            dq.push_back(arg2);
                    }
                    // Such things as a+0 or 0+a
                    else if(it->str() == "+" && ((arg2.is_constant() && arg2.eval() == static_cast<T>(0.0)) ||
                                                 (arg1.is_constant() && arg1.eval() == static_cast<T>(0.0) &&
                                                  !arg2.is_operator() && !arg2.is_function())))
                    {
                        dq.pop_back();
                        if(arg2.is_constant() && arg2.eval() == static_cast<T>(0.0))
                            dq.push_back(arg1);
                        else
                            dq.push_back(arg2);
                    }
                    // Such things as a-0
                    else if(it->str() == "-" && arg2.is_constant() && arg2.eval() == static_cast<T>(0.0))
                    {
                        dq.pop_back();
                        dq.push_back(arg1);
                    }
                    // Nothing...
                    else
                    {
                        dq.push_back(arg2);
                        dq.push_back(*it);
                    }
                }
                else
                {
                    dq.push_back(*it);
                }
            }

            if(expression_objects.size() > dq.size())
            {
                expression_objects.clear();
                expression_objects.reserve(dq.size());
                while(!dq.empty())
                {
                    expression_objects.push_back(dq.front());
                    dq.pop_front();
                }
                was_changed = true;
            }
            else
            {
                dq.clear();
            }
        }
        while(was_changed);
        return true;
    }

    bool calculate(T & result)
    {
        using namespace std;
        using namespace parser_internal;

        if(!is_parsed())
        {
            error_string = "Not parsed!";
            return false;
        }

        stack<T> st;

        for(typename vector<parser_object<T> >::const_iterator it = expression_objects.begin(); it != expression_objects.end(); ++it)
        {
            if(it->is_constant())
            {
                st.push(it->eval());
            }
            else if(it->is_variable())
            {
                typename map<string, T>::const_iterator itc = constants.find(it->str());
                if(itc != constants.end())
                {
                    st.push(itc->second);
                }
                else
                {
                    stringstream sst;
                    sst << "Constant `" << it->str() << "` must be defined!";
                    error_string = sst.str();
                    return false;
                }
            }
            else if(it->is_operator())
            {
                T arg2 = st.top();
                st.pop();
                T arg1 = st.top();
                st.pop();
                T val = it->eval(arg1, arg2);
                st.push(val);
            }
            else if(it->is_function())
            {
                T arg = st.top();
                st.pop();
                T val = it->eval(arg);
                st.push(val);
            }
        }

        if(st.size() != 1)
        {
            stringstream sst;
            sst << "Stack size equal " << st.size();
            error_string = sst.str();
            return false;
        }
        result = st.top();
        return true;
    }

    void debug_print() const
    {
        using namespace std;
        using namespace parser_internal;
        if(expression.size() > expression_objects.size())
            for(vector<string>::const_iterator it = expression.begin(); it != expression.end(); ++it)
                cout << * it << ' ';
        else
            for(typename vector<parser_object<T> >::const_iterator it = expression_objects.begin(); it != expression_objects.end(); ++it)
                cout << it->str() << ' ';
        cout << endl;
    }
};

#endif // PARSER_H
