#include "config.h"
#include <fstream>
#include <cassert>

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
    string result;
    for(string::const_iterator it = str.begin(); it != str.end(); ++it)
    {
        char c = * it;
        if(c >= 'A' && c <= 'Z')
            c -= 'A' - 'a';
        result += c;
    }
    return result;
}

// ============================================================================

evaluator::evaluator()
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
    filename_phys = "phys.txt";
    filename_slae = "";

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
                            else cerr << "[Config] Unsupported param \"" << param<< "\" in section \"" << section
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
                array_t<parser<complex<double> >, 3> * curr_parser;
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
                            if(param == "x")
                            {
                                if(!(* curr_parser)[0].parse(value) || !(* curr_parser)[0].simplify())
                                {
                                    cerr << "Error in " << __FILE__ << ":" << __LINE__
                                         << " while parsing " << value << endl;
                                    return init(false);
                                }
                            }
                            else if(param == "y")
                            {
                                if(!(* curr_parser)[1].parse(value) || !(* curr_parser)[1].simplify())
                                {
                                    cerr << "Error in " << __FILE__ << ":" << __LINE__
                                         << " while parsing " << value << endl;
                                    return init(false);
                                }
                            }
                            else if(param == "z")
                            {
                                if(!(* curr_parser)[2].parse(value) || !(* curr_parser)[2].simplify())
                                {
                                    cerr << "Error in " << __FILE__ << ":" << __LINE__
                                         << " while parsing " << value << endl;
                                    return init(false);
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
                            else cerr << "[Config] Unsupported param \"" << param<< "\" in section \"" << section
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
                            else cerr << "[Config] Unsupported param \"" << param<< "\" in section \"" << section
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
