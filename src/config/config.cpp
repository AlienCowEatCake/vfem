#include "config.h"
#include <fstream>
#include <cassert>
#include <cctype>
#include <algorithm>

// ============================================================================

string trim(const string & str)
{
    size_t start = str.find_first_not_of(" \t\f\v\n\r");
    size_t stop = str.find_last_not_of(" \t\f\v\n\r");
    if(start == string::npos || stop == string::npos)
        return "";
    return str.substr(start, stop - start + 1);
}

string to_lowercase(const string & str)
{
    string result = str;
    transform(str.begin(), str.end(), result.begin(), ::tolower);
    return result;
}

// ============================================================================

evaluator3::evaluator3()
{
    for(size_t i = 0; i < 3; i++)
        default_value[i].parse("0.0");
}

// ============================================================================

config_type::config_type()
{
    load_defaults();
}

// Загрузка значений по-умолчанию
void config_type::load_defaults()
{
    basis.order = 1;
    basis.type = 2;
    eps_slae = 1e-10;
    eps_slae_bound = 1e-14;
    gamma_v_cycle_0 = 0.01;
    gamma_v_cycle_full = 0.05;
    gamma_v_cycle_ker = 0.01;
    max_iter_v_cycle_local = 500;
    filename_mesh = "mesh.msh";
    filename_phys = "phys.ini";
    filename_pml = "config_pml.ini";
    filename_slae = "";
    jit_type = evaluator3::JIT_DISABLE;

    analytical_enabled = false;
    boundary_enabled = false;
    right_enabled = false;

    init(true);
}

// Пост-загрузочная инициализация
bool config_type::init(bool status)
{
    assert(basis.order >= 1 && basis.order <= 2);
    assert(basis.type == 1 || basis.type == 2);
    // Базис первого неполного порядка
    if(basis.order == 1 && basis.type == 1)
    {
        basis.tet_bf_num = 6;
        basis.tet_ker_bf_num = 4;
        basis.tr_bf_num = 3;
        basis.tr_ker_bf_num = 3;
    }
    // Базис первого полного порядка
    else if(basis.order == 1 && basis.type == 2)
    {
        basis.tet_bf_num = 12;
        basis.tet_ker_bf_num = 10;
        basis.tr_bf_num = 6;
        basis.tr_ker_bf_num = 6;
    }
    // Базис второго неполного порядка
    else if(basis.order == 2 && basis.type == 1)
    {
        basis.tet_bf_num = 20;
        basis.tet_ker_bf_num = 10;
        basis.tr_bf_num = 8;
        basis.tr_ker_bf_num = 6;
    }
    // Базис второго полного порядка
    else if(basis.order == 2 && basis.type == 2)
    {
        basis.tet_bf_num = 30;
        basis.tet_ker_bf_num = 20;
        basis.tr_bf_num = 12;
        basis.tr_ker_bf_num = 10;
    }
    return status;
}

// Загрузка значений из файла
bool config_type::load(const string & filename)
{
    cout << "Reading config file ..." << endl;

    ifstream ifs(filename.c_str());
    if(!ifs.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return init(false);
    }

    string line;
    getline(ifs, line);
    line = trim(line);
    do
    {
        if(line.length() > 1 && line[0] == '[')
        {
            line = to_lowercase(line.substr(1, line.length() - 2));
            string section = line, subsection;
            size_t dot_pos = line.find_first_of(".");
            if(dot_pos != string::npos)
            {
                section = line.substr(0, dot_pos);
                subsection = line.substr(dot_pos + 1);
            }
            //cout << "section: " << section << "\nsubsection: " << subsection << endl;

            if(section == "vfem")
            {
                do
                {
                    getline(ifs, line);
                    line = trim(line);
                    if(line.length() > 1 && line[0] != ';')
                    {
                        size_t eq_pos = line.find_first_of("=");
                        if(eq_pos != string::npos)
                        {
                            string param = to_lowercase(trim(line.substr(0, eq_pos)));
                            string value = trim(line.substr(eq_pos + 1));
                            if(value.length() > 1 && value[0] == '\"')
                                value = trim(value.substr(1, value.length() - 2));
                            stringstream sst(value);
                            if(param == "basis_order")                  sst >> basis.order;
                            else if(param == "basis_type")              sst >> basis.type;
                            else if(param == "eps_slae")                sst >> eps_slae;
                            else if(param == "eps_slae_bound")          sst >> eps_slae_bound;
                            else if(param == "gamma_v_cycle_full")      sst >> gamma_v_cycle_full;
                            else if(param == "gamma_v_cycle_ker")       sst >> gamma_v_cycle_ker;
                            else if(param == "gamma_v_cycle_0")         sst >> gamma_v_cycle_0;
                            else if(param == "max_iter_v_cycle_local")  sst >> max_iter_v_cycle_local;
                            else if(param == "filename_mesh")           sst >> filename_mesh;
                            else if(param == "filename_phys")           sst >> filename_phys;
                            else if(param == "filename_slae")           sst >> filename_slae;
                            else if(param == "filename_pml")            sst >> filename_pml;
                            else if(param == "jit_type")
                            {
                                value = to_lowercase(value);
                                if     (value == "inline")  jit_type = evaluator3::JIT_INLINE;
                                else if(value == "extcall") jit_type = evaluator3::JIT_EXTCALL;
                                else if(value == "disable") jit_type = evaluator3::JIT_DISABLE;
                                else cerr << "[Config] Unknown JIT-compiler type \"" << value << "\" in section \""
                                          << section << (subsection == "" ? string("") : (string(".") + subsection))
                                          << "\"" << endl;
                            }
                            else cerr << "[Config] Unsupported param \"" << param << "\" in section \"" << section
                                      << (subsection == "" ? string("") : (string(".") + subsection)) << "\"" << endl;
                            //cout << "  param = " << param << endl;
                            //cout << "  value = " << value << endl;
                        }
                    }
                }
                while(ifs.good() && !(line.length() > 1 && line[0] == '['));
            }

            else if(section == "boundary" || section == "right" || section == "analytical")
            {
                array_t<evaluator<complex<double> >, 3> * curr_parser = NULL;
                if(section == "boundary")
                {
                    if(subsection == "")
                        curr_parser = & boundary.default_value;
                    else
                    {
                        stringstream sst(subsection);
                        size_t bound_phys;
                        sst >> bound_phys;
                        curr_parser = & boundary.values[bound_phys];
                    }
                }
                else if(section == "right")
                {
                    if(subsection == "")
                        curr_parser = & right.default_value;
                    else
                    {
                        stringstream sst(subsection);
                        size_t bound_phys;
                        sst >> bound_phys;
                        curr_parser = & right.values[bound_phys];
                    }
                }
                else if(section == "analytical")
                {
                    if(subsection == "")
                        curr_parser = & analytical.default_value;
                    else
                    {
                        stringstream sst(subsection);
                        size_t bound_phys;
                        sst >> bound_phys;
                        curr_parser = & analytical.values[bound_phys];
                    }
                }

                do
                {
                    getline(ifs, line);
                    line = trim(line);
                    if(line.length() > 1 && line[0] != ';')
                    {
                        size_t eq_pos = line.find_first_of("=");
                        if(eq_pos != string::npos)
                        {
                            string param = to_lowercase(trim(line.substr(0, eq_pos)));
                            string value = trim(line.substr(eq_pos + 1));
                            if(value.length() > 1 && value[0] == '\"')
                                value = trim(value.substr(1, value.length() - 2));
                            if(param == "x" || param == "y" || param == "z")
                            {
                                evaluator<complex<double> > * curr_parser_var = NULL;
                                if(param == "x")      curr_parser_var = & ((* curr_parser)[0]);
                                else if(param == "y") curr_parser_var = & ((* curr_parser)[1]);
                                else if(param == "z") curr_parser_var = & ((* curr_parser)[2]);
                                if(!curr_parser_var->parse(value) || !curr_parser_var->simplify())
                                {
                                    cerr << "Error in " << __FILE__ << ":" << __LINE__
                                         << " while parsing " << value << " ("
                                         << curr_parser_var->get_error() << ")" << endl;
                                    return init(false);
                                }
                                switch(jit_type)
                                {
                                case evaluator3::JIT_EXTCALL:
                                    if(!curr_parser_var->compile_extcall())
                                        cerr << "Error in " << __FILE__ << ":" << __LINE__
                                             << " while compiling (extcall) " << value << " ("
                                             << curr_parser_var->get_error() << ")" << endl;
                                    break;
                                case evaluator3::JIT_INLINE:
                                    if(!curr_parser_var->compile_inline())
                                        cerr << "Error in " << __FILE__ << ":" << __LINE__
                                             << " while compiling (inline) " << value << " ("
                                             << curr_parser_var->get_error() << ")" << endl;
                                    break;
                                default:
                                    break;
                                }
                            }
                            else if(param == "enabled")
                            {
                                if(section == "analytical")
                                {
                                    value = to_lowercase(value);
                                    if(value == "yes" || value == "true" || value == "1")
                                        analytical_enabled = true;
                                    else
                                        analytical_enabled = false;
                                }
                                else if(section == "boundary")
                                {
                                    value = to_lowercase(value);
                                    if(value == "yes" || value == "true" || value == "1")
                                        boundary_enabled = true;
                                    else
                                        boundary_enabled = false;
                                }
                                else if(section == "right")
                                {
                                    value = to_lowercase(value);
                                    if(value == "yes" || value == "true" || value == "1")
                                        right_enabled = true;
                                    else
                                        right_enabled = false;
                                }
                            }
                            else cerr << "[Config] Unsupported param \"" << param << "\" in section \"" << section
                                      << (subsection == "" ? string("") : (string(".") + subsection)) << "\"" << endl;
                            //cout << "  param = " << param << endl;
                            //cout << "  value = " << value << endl;
                        }
                    }
                }
                while(ifs.good() && !(line.length() > 1 && line[0] == '['));
            }

            else if(section == "postprocessing")
            {
                stringstream sst(subsection);
                size_t index;
                sst >> index;
                postprocessor * p = &(post[index]);

                do
                {
                    getline(ifs, line);
                    line = trim(line);
                    if(line.length() > 1 && line[0] != ';')
                    {
                        size_t eq_pos = line.find_first_of("=");
                        if(eq_pos != string::npos)
                        {
                            string param = to_lowercase(trim(line.substr(0, eq_pos)));
                            string value = trim(line.substr(eq_pos + 1));
                            if(value.length() > 1 && value[0] == '\"')
                                value = trim(value.substr(1, value.length() - 2));
                            stringstream sst(value);
                            if(param == "filename")             sst >> p->filename;
                            else if(param == "slice_var_name")  sst >> p->param_2d.slice_var_name;
                            else if(param == "slice_var_value") sst >> p->param_2d.slice_var_value;
                            else if(param == "var1_name")       sst >> p->param_2d.var1_name;
                            else if(param == "var1_from")       sst >> p->param_2d.var1_from;
                            else if(param == "var1_to")         sst >> p->param_2d.var1_to;
                            else if(param == "var1_num")        sst >> p->param_2d.var1_num;
                            else if(param == "var2_name")       sst >> p->param_2d.var2_name;
                            else if(param == "var2_from")       sst >> p->param_2d.var2_from;
                            else if(param == "var2_to")         sst >> p->param_2d.var2_to;
                            else if(param == "var2_num")        sst >> p->param_2d.var2_num;
                            else if(param == "line_var1_name")  sst >> p->param_1d.line_var1_name;
                            else if(param == "line_var1_value") sst >> p->param_1d.line_var1_value;
                            else if(param == "line_var2_name")  sst >> p->param_1d.line_var2_name;
                            else if(param == "line_var2_value") sst >> p->param_1d.line_var2_value;
                            else if(param == "var_name")        sst >> p->param_1d.var_name;
                            else if(param == "var_from")        sst >> p->param_1d.var_from;
                            else if(param == "var_to")          sst >> p->param_1d.var_to;
                            else if(param == "var_num")         sst >> p->param_1d.var_num;
                            else if(param == "timestamp")
                            {
                                value = to_lowercase(value);
                                if(value == "yes" || value == "true" || value == "1")
                                    p->timestamp = true;
                                else
                                    p->timestamp = false;
                            }
                            else if(param == "type")
                            {
                                if(value == "1d") p->type = 1;
                                if(value == "2d") p->type = 2;
                                if(value == "3d") p->type = 3;
                            }
                            else cerr << "[Config] Unsupported param \"" << param << "\" in section \"" << section
                                      << (subsection == "" ? string("") : (string(".") + subsection)) << "\"" << endl;
                            //cout << "  param = " << param << endl;
                            //cout << "  value = " << value << endl;
                        }
                    }
                }
                while(ifs.good() && !(line.length() > 1 && line[0] == '['));
            }

        }
        else
        {
            getline(ifs, line);
            line = trim(line);
        }
    }
    while(ifs.good());

    ifs.close();
    return init(true);
}

bool config_type::load_pml(const string & filename)
{
    cout << "Reading PML config file ..." << endl;

    ifstream ifs(filename.c_str());
    if(!ifs.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return false;
    }

    // Умолчательные значения зададим пока огромными числами
    complex<double> chi_def(DBL_MAX, DBL_MAX);
    double m_def = DBL_MAX;
    double width_def = DBL_MAX;
    phys_pml.x0 = DBL_MAX;
    phys_pml.x1 = -DBL_MAX;
    phys_pml.y0 = DBL_MAX;
    phys_pml.y1 = -DBL_MAX;
    phys_pml.z0 = DBL_MAX;
    phys_pml.z1 = -DBL_MAX;

    string line;
    getline(ifs, line);
    line = trim(line);
    do
    {
        if(line.length() > 1 && line[0] == '[')
        {
            line = to_lowercase(line.substr(1, line.length() - 2));
            string section = line, subsection;
            size_t dot_pos = line.find_first_of(".");
            if(dot_pos != string::npos)
            {
                section = line.substr(0, dot_pos);
                subsection = line.substr(dot_pos + 1);
            }
            //cout << "section: " << section << "\nsubsection: " << subsection << endl;

            if(section == "pml")
            {
                if(subsection == "")
                {
                    do
                    {
                        getline(ifs, line);
                        line = trim(line);
                        if(line.length() > 1 && line[0] != ';')
                        {
                            size_t eq_pos = line.find_first_of("=");
                            if(eq_pos != string::npos)
                            {
                                string param = to_lowercase(trim(line.substr(0, eq_pos)));
                                string value = trim(line.substr(eq_pos + 1));
                                if(value.length() > 1 && value[0] == '\"')
                                    value = trim(value.substr(1, value.length() - 2));
                                stringstream sst(value);
                                double tmp;
                                sst >> tmp;
                                if(param == "chi_real")         chi_def.real(tmp);
                                else if(param == "chi_imag")    chi_def.imag(tmp);
                                else if(param == "m")           m_def = tmp;
                                else if(param == "width")       width_def = tmp;
                                else if(param == "x0")          phys_pml.x0 = tmp;
                                else if(param == "x1")          phys_pml.x1 = tmp;
                                else if(param == "y0")          phys_pml.y0 = tmp;
                                else if(param == "y1")          phys_pml.y1 = tmp;
                                else if(param == "z0")          phys_pml.z0 = tmp;
                                else if(param == "z1")          phys_pml.z1 = tmp;
                                else cerr << "[PML Config] Unsupported param \"" << param << "\" in section \"" << section
                                          << (subsection == "" ? string("") : (string(".") + subsection)) << "\"" << endl;
                                //cout << "  param = " << param << endl;
                                //cout << "  value = " << value << endl;
                            }
                        }
                    }
                    while(ifs.good() && !(line.length() > 1 && line[0] == '['));
                }
                else
                {
                    stringstream sst(subsection);
                    size_t pml_phys;
                    sst >> pml_phys;
                    pml_config_parameter * curr = NULL;
                    if(phys_pml.params.find(pml_phys) != phys_pml.params.end())
                    {
                        curr = & (phys_pml.params[pml_phys]);
                    }
                    else
                    {
                        curr = & (phys_pml.params[pml_phys]);
                        curr->chi = chi_def;
                        curr->m = m_def;
                        curr->width = width_def;
                    }

                    do
                    {
                        getline(ifs, line);
                        line = trim(line);
                        if(line.length() > 1 && line[0] != ';')
                        {
                            size_t eq_pos = line.find_first_of("=");
                            if(eq_pos != string::npos)
                            {
                                string param = to_lowercase(trim(line.substr(0, eq_pos)));
                                string value = trim(line.substr(eq_pos + 1));
                                if(value.length() > 1 && value[0] == '\"')
                                    value = trim(value.substr(1, value.length() - 2));
                                stringstream sst(value);
                                double tmp;
                                sst >> tmp;
                                if(param == "chi_real")         curr->chi.real(tmp);
                                else if(param == "chi_imag")    curr->chi.imag(tmp);
                                else if(param == "m")           curr->m = tmp;
                                else if(param == "width")       curr->width = tmp;
                                else cerr << "[PML Config] Unsupported param \"" << param << "\" in section \"" << section
                                          << (subsection == "" ? string("") : (string(".") + subsection)) << "\"" << endl;
                                //cout << "  param = " << param << endl;
                                //cout << "  value = " << value << endl;
                            }
                        }
                    }
                    while(ifs.good() && !(line.length() > 1 && line[0] == '['));
                }
            }
        }
        else
        {
            getline(ifs, line);
            line = trim(line);
        }
    }
    while(ifs.good());

    ifs.close();

    // Проверим введенное на адекватность
    pml_config_parameter def;
    double big_num = DBL_MAX * 0.999;
    if(chi_def.real() < big_num) def.chi.real(chi_def.real());
    if(chi_def.imag() < big_num) def.chi.imag(chi_def.real());
    if(m_def < big_num)          def.m = m_def;
    if(width_def < big_num)      def.width = width_def;
    for(map<size_t, pml_config_parameter>::iterator it = phys_pml.params.begin(); it != phys_pml.params.end(); ++it)
    {
        if(it->second.chi.real() >= big_num) it->second.chi.real(def.chi.real());
        if(it->second.chi.imag() >= big_num) it->second.chi.imag(def.chi.imag());
        if(it->second.m >= big_num)          it->second.m = def.m;
        if(it->second.width >= big_num)      it->second.width = def.width;
    }

    return true;
}
