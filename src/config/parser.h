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

template<typename T>
class parser
{
protected:
    std::vector<std::string> expression;
    std::set<std::string> functions;
    std::map<std::string, T> constants;
    std::map<char, unsigned short int> operators;
    bool status;
    std::string error_string;

    bool is_complex;
    T imag_value;

    void init()
    {
        using namespace std;

        if(typeid(T) == typeid(complex<float>) ||
           typeid(T) == typeid(complex<double>) ||
           typeid(T) == typeid(complex<long double>))
        {
            is_complex = true;
            imag_value = sqrt(static_cast<T>(-1.0));
            functions.insert("imag");
            functions.insert("real");
        }
        else
        {
            is_complex = false;
        }

        functions.insert("sin");
        functions.insert("cos");
        functions.insert("tan");
        functions.insert("asin");
        functions.insert("acos");
        functions.insert("atan");
        functions.insert("sinh");
        functions.insert("cosh");
        functions.insert("tanh");
        functions.insert("asinh");
        functions.insert("acosh");
        functions.insert("atanh");
        functions.insert("log");
        functions.insert("log10");
        functions.insert("abs");
        functions.insert("exp");
        functions.insert("sqrt");

        operators['+'] = 1;
        operators['-'] = 1;
        operators['*'] = 2;
        operators['/'] = 2;
        operators['^'] = 3;
    }

    void init_const()
    {
        constants["pi"] = static_cast<T>(3.14159265358979323846264338327950);
        constants["e"]  = static_cast<T>(2.71828182845904523536028747135266);
        if(is_complex)
        {
            constants["i"] = imag_value;
            constants["j"] = imag_value;
        }
    }

    bool calc_operator(const std::string & oper, const T & larg, const T & rarg, T & result)
    {
        using namespace std;
        switch(oper[0])
        {
        case '+':
            result = larg + rarg;
            return true;
        case '-':
            result = larg - rarg;
            return true;
        case '*':
            result = larg * rarg;
            return true;
        case '/':
            result = larg / rarg;
            return true;
        case '^':
            result = pow(larg, rarg);
            return true;
        }
        stringstream sst;
        sst << "Unknown operator " << oper;
        error_string = sst.str();
        return false;
    }

    template<typename U>
    bool calc_function(const std::string & func, const std::complex<U> & arg, T & result)
    {
        using namespace std;
        if(func == "sin")
        {
            result = sin(arg);
            return true;
        }
        if(func == "cos")
        {
            result = cos(arg);
            return true;
        }
        if(func == "tan")
        {
            result = tan(arg);
            return true;
        }
        if(func == "sinh")
        {
            result = sinh(arg);
            return true;
        }
        if(func == "cosh")
        {
            result = cosh(arg);
            return true;
        }
        if(func == "tanh")
        {
            result = tanh(arg);
            return true;
        }
        if(func == "log")
        {
            result = log(arg);
            return true;
        }
        if(func == "log10")
        {
            result = log10(arg);
            return true;
        }
        if(func == "abs")
        {
            result = abs(arg);
            return true;
        }
        if(func == "sqrt")
        {
            result = sqrt(arg);
            return true;
        }
        if(func == "exp")
        {
            result = exp(arg);
            return true;
        }
        if(func == "imag")
        {
            result = arg.imag();
            return true;
        }
        if(func == "real")
        {
            result = arg.real();
            return true;
        }
        if(fabs(arg.imag()) < numeric_limits<U>::epsilon())
        {
            if(func == "asin")
            {
                result = asin(arg.real());
                return true;
            }
            if(func == "acos")
            {
                result = acos(arg.real());
                return true;
            }
            if(func == "atan")
            {
                result = atan(arg.real());
                return true;
            }
            if(func == "asinh")
            {
                result = asinh(arg.real());
                return true;
            }
            if(func == "acosh")
            {
                result = acosh(arg.real());
                return true;
            }
            if(func == "atanh")
            {
                result = atanh(arg.real());
                return true;
            }
        }
        stringstream sst;
        sst << "Unknown function " << func;
        error_string = sst.str();
        return false;
    }

    template<typename U>
    bool calc_function(const std::string & func, const U & arg, T & result)
    {
        using namespace std;
        if(func == "sin")
        {
            result = sin(arg);
            return true;
        }
        if(func == "cos")
        {
            result = cos(arg);
            return true;
        }
        if(func == "tan")
        {
            result = tan(arg);
            return true;
        }
        if(func == "asin")
        {
            result = asin(arg);
            return true;
        }
        if(func == "acos")
        {
            result = acos(arg);
            return true;
        }
        if(func == "atan")
        {
            result = atan(arg);
            return true;
        }
        if(func == "sinh")
        {
            result = sinh(arg);
            return true;
        }
        if(func == "cosh")
        {
            result = cosh(arg);
            return true;
        }
        if(func == "tanh")
        {
            result = tanh(arg);
            return true;
        }
        if(func == "asinh")
        {
            result = asinh(arg);
            return true;
        }
        if(func == "acosh")
        {
            result = acosh(arg);
            return true;
        }
        if(func == "atanh")
        {
            result = atanh(arg);
            return true;
        }
        if(func == "log")
        {
            result = log(arg);
            return true;
        }
        if(func == "log10")
        {
            result = log10(arg);
            return true;
        }
        if(func == "abs")
        {
            result = fabs(arg);
            return true;
        }
        if(func == "sqrt")
        {
            result = sqrt(arg);
            return true;
        }
        if(func == "exp")
        {
            result = exp(arg);
            return true;
        }
        stringstream sst;
        sst << "Unknown function " << func;
        error_string = sst.str();
        return false;
    }

public:
    parser()
    {
        init();
        init_const();
        status = false;
    }

    parser(const std::string & str)
    {
        init();
        init_const();
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
        constants.clear();
        init_const();
    }

    bool is_parsed() const
    {
        return status;
    }

    bool parse(const std::string & str)
    {
        using namespace std;

        expression.clear();
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
                    sym = *it;
                }
                if(it != str.end() && (sym == '.' || sym == ','))
                {
                    a.push_back('.');
                    ++it;
                    sym = *it;
                    while(it != str.end() && sym >= '0' && sym <= '9')
                    {
                        a.push_back(sym);
                        ++it;
                        sym = *it;
                    }
                }
                if(it != str.end() && (sym == 'e' || sym == 'E' || sym == 'd' || sym == 'D'))
                {
                    a.push_back('e');
                    ++it;
                    sym = *it;
                    if(it != str.end() && (sym == '-' || sym == '+'))
                    {
                        a.push_back(sym);
                        ++it;
                        sym = *it;
                    }
                    while(it != str.end() && sym >= '0' && sym <= '9')
                    {
                        a.push_back(sym);
                        ++it;
                        sym = *it;
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
                          operators[sym] <= operators[op])
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

        return status;
    }

    bool simplify()
    {
        using namespace std;

        if(!is_parsed())
        {
            error_string = "Not parsed!";
            return false;
        }

        bool was_changed;
        do
        {
            deque<string> dq;
            was_changed = false;

            for(vector<string>::iterator it = expression.begin(); it != expression.end(); ++it)
            {
                if(constants.find(*it) == constants.end() && operators.find((*it)[0]) != operators.end())
                {
                    string arg2 = dq.back();
                    if(constants.find(arg2) != constants.end())
                    {
                        dq.pop_back();
                        string arg1 = dq.back();
                        if(constants.find(arg1) != constants.end())
                        {
                            dq.pop_back();
                            T val;
                            if(!calc_operator(*it, constants[arg1], constants[arg2], val))
                                return false;
                            stringstream sst;
                            sst.precision(17);
                            sst.setf(ios::scientific);
                            sst << val;
                            string sst_st = sst.str();
                            constants[sst_st] = val;
                            dq.push_back(sst_st);
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
                else if(functions.find(*it) != functions.end())
                {
                    string arg = dq.back();
                    if(constants.find(arg) != constants.end())
                    {
                        dq.pop_back();
                        T val;
                        if(!calc_function(*it, constants[arg], val))
                            return false;
                        stringstream sst;
                        sst.precision(17);
                        sst.setf(ios::scientific);
                        sst << val;
                        string sst_st = sst.str();
                        constants[sst_st] = val;
                        dq.push_back(sst_st);
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

            if(expression.size() > dq.size())
            {
                expression.clear();
                expression.reserve(dq.size());
                while(!dq.empty())
                {
                    expression.push_back(dq.front());
                    dq.pop_front();
                }
                was_changed = true;
            }
            else
            {
                dq.clear();
            }

            for(vector<string>::iterator it = expression.begin(); it != expression.end(); ++it)
            {
                if(constants.find(*it) == constants.end() && operators.find((*it)[0]) != operators.end())
                {
                    string arg2 = dq.back();
                    dq.pop_back();
                    string arg1 = dq.back();
                    // Such things as a*0 or 0*a
                    if(*it == "*" && ((constants.find(arg2) != constants.end() && constants[arg2] == static_cast<T>(0.0)) ||
                                      (constants.find(arg1) != constants.end() && constants[arg1] == static_cast<T>(0.0) &&
                                       operators.find(arg2[0]) == operators.end() &&
                                       functions.find(arg2) == functions.end())))
                    {
                        dq.pop_back();
                        if(constants.find(arg2) != constants.end() && constants[arg2] == static_cast<T>(0.0))
                            dq.push_back(arg2);
                        else
                            dq.push_back(arg1);
                    }
                    // Such things as a*1 or 1*a
                    else if(*it == "*" && ((constants.find(arg2) != constants.end() && constants[arg2] == static_cast<T>(1.0)) ||
                                           (constants.find(arg1) != constants.end() && constants[arg1] == static_cast<T>(1.0) &&
                                            operators.find(arg2[0]) == operators.end() &&
                                            functions.find(arg2) == functions.end())))
                    {
                        dq.pop_back();
                        if(constants.find(arg2) != constants.end() && constants[arg2] == static_cast<T>(1.0))
                            dq.push_back(arg1);
                        else
                            dq.push_back(arg2);
                    }
                    // Such things as a+0 or 0+a
                    else if(*it == "+" && ((constants.find(arg2) != constants.end() && constants[arg2] == static_cast<T>(0.0)) ||
                                           (constants.find(arg1) != constants.end() && constants[arg1] == static_cast<T>(0.0) &&
                                            operators.find(arg2[0]) == operators.end() &&
                                            functions.find(arg2) == functions.end())))
                    {
                        dq.pop_back();
                        if(constants.find(arg2) != constants.end() && constants[arg2] == static_cast<T>(0.0))
                            dq.push_back(arg1);
                        else
                            dq.push_back(arg2);
                    }
                    // Such things as a-0
                    else if(*it == "-" && constants.find(arg2) != constants.end() && constants[arg2] == static_cast<T>(0.0))
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

            if(expression.size() > dq.size())
            {
                expression.clear();
                expression.reserve(dq.size());
                while(!dq.empty())
                {
                    expression.push_back(dq.front());
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

        if(!is_parsed())
        {
            error_string = "Not parsed!";
            return false;
        }

        stack<T> st;

        for(vector<string>::const_iterator it = expression.begin(); it != expression.end(); ++it)
        {
            if(constants.find(*it) != constants.end())
            {
                st.push((*constants.find(*it)).second);
            }
            else if(operators.find((*it)[0]) != operators.end())
            {
                T arg2 = st.top();
                st.pop();
                T arg1 = st.top();
                st.pop();
                T val;
                if(!calc_operator(*it, arg1, arg2, val))
                    return false;
                st.push(val);
            }
            else if(functions.find(*it) != functions.end())
            {
                T arg = st.top();
                st.pop();
                T val;
                if(!calc_function(*it, arg, val))
                    return false;
                st.push(val);
            }
            else
            {
                error_string = "Constants must be defined!";
                return false;
            }
        }

        if(st.size() != 1)
        {
            error_string = "Internal error!";
            cerr << "[parser] Internal error: Stack size equal " << st.size() << endl;
            return false;
        }
        result = st.top();
        return true;
    }

    void debug_print() const
    {
        using namespace std;
        for(vector<string>::const_iterator it = expression.begin(); it != expression.end(); ++it)
            cout << * it << ' ';
        cout << endl;
    }
};

#endif // PARSER_H
