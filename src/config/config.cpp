#include "config.h"

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
    jit_types_table[evaluator3::JIT_DISABLE] = "disable";
    jit_types_table[evaluator3::JIT_INLINE]  = "inline";
    jit_types_table[evaluator3::JIT_EXTCALL] = "extcall";
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
    v_cycle_enabled = false;
    max_iter = 15000;
    max_iter_v_cycle_local = 500;
    filename_mesh = "mesh.msh";
    filename_phys = "phys.ini";
    filename_pml = "";
    filename_slae = "";
    jit_type = evaluator3::JIT_DISABLE;

    analytical_enabled = false;
    boundary_enabled = false;
    right_enabled = false;

    // Если это не исправлено при чтении конфига с PML,
    // значит пусть вся область считается как область без PML
    phys_pml.x0 = -DBL_MAX;
    phys_pml.x1 = DBL_MAX;
    phys_pml.y0 = -DBL_MAX;
    phys_pml.y1 = DBL_MAX;
    phys_pml.z0 = -DBL_MAX;
    phys_pml.z1 = DBL_MAX;
    phys_pml.params.clear();

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
    // Если точность решения СЛАУ больше точности первичного
    // уточнения решения, значит следует выключить V-цикл
    if(eps_slae >= gamma_v_cycle_0)
    {
        eps_slae = gamma_v_cycle_0;
        v_cycle_enabled = false;
    }
    return status;
}

// Загрузка значений из файла
bool config_type::load(const string & filename)
{
    cout << "Reading config file ..." << endl;

    inifile cfg_file(filename);
    if(!cfg_file.good())
    {
        cout << "[Config] Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return init(false);
    }

    // Секция VFEM

    basis.order             = cfg_file.get("vfem", "", "basis_order",            basis.order);
    basis.type              = cfg_file.get("vfem", "", "basis_type",             basis.type);
    eps_slae                = cfg_file.get("vfem", "", "eps_slae",               eps_slae);
    eps_slae_bound          = cfg_file.get("vfem", "", "eps_slae_bound",         eps_slae_bound);
    gamma_v_cycle_full      = cfg_file.get("vfem", "", "gamma_v_cycle_full",     gamma_v_cycle_full);
    gamma_v_cycle_ker       = cfg_file.get("vfem", "", "gamma_v_cycle_ker",      gamma_v_cycle_ker);
    gamma_v_cycle_0         = cfg_file.get("vfem", "", "gamma_v_cycle_0",        gamma_v_cycle_0);
    max_iter_v_cycle_local  = cfg_file.get("vfem", "", "max_iter_v_cycle_local", max_iter_v_cycle_local);
    max_iter                = cfg_file.get("vfem", "", "max_iter",               max_iter);
    filename_mesh           = cfg_file.get("vfem", "", "filename_mesh",          filename_mesh);
    filename_phys           = cfg_file.get("vfem", "", "filename_phys",          filename_phys);
    filename_slae           = cfg_file.get("vfem", "", "filename_slae",          filename_slae);
    filename_pml            = cfg_file.get("vfem", "", "filename_pml",           filename_pml);
    string jit_type_str     = cfg_file.get("vfem", "", "jit_type",               jit_types_table[jit_type]);
    v_cycle_enabled         = cfg_file.get("vfem", "", "v_cycle_enabled",        v_cycle_enabled);

    int tet_order           = cfg_file.get("vfem", "", "tet_integration_order",  (int)-1);
    if(tet_order > 0)
        basis.tet_int.init((size_t)tet_order);

    int tr_order            = cfg_file.get("vfem", "", "tr_integration_order",   (int)-1);
    if(tr_order > 0)
        basis.tr_int.init((size_t)tr_order);

    jit_type_str = to_lowercase(jit_type_str);
    for(map<evaluator3::jit_types, string>::iterator it = jit_types_table.begin(); it != jit_types_table.end(); ++it)
    {
        if(it->second == jit_type_str)
            jit_type = it->first;
    }

    // Секция Boundary

    list<size_t> boundaries = cfg_file.enumerate("boundary", (size_t)1);
    for(list<size_t>::iterator it = boundaries.begin(); it != boundaries.end(); ++it)
    {
        size_t subsection = * it;
        array_t<evaluator<complex<double> >, 3> * curr_parser = NULL;
        if(subsection == 0)
            curr_parser = & boundary.default_value;
        else
            curr_parser = & boundary.values[subsection];

        curr_parser[0] = cfg_file.get("boundary", subsection, )
    }






            else if(section == "boundary" || section == "right" || section == "analytical")
            {

                if(section == "boundary")
                {
                    if(subsection.empty())
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
                    if(subsection.empty())
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
                    if(subsection.empty())
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
                                    cout << "Error in " << __FILE__ << ":" << __LINE__
                                         << " while parsing " << value << " ("
                                         << curr_parser_var->get_error() << ")" << endl;
                                    return init(false);
                                }
                                switch(jit_type)
                                {
                                case evaluator3::JIT_EXTCALL:
                                    if(!curr_parser_var->compile_extcall())
                                        cout << "Error in " << __FILE__ << ":" << __LINE__
                                             << " while compiling (extcall) " << value << " ("
                                             << curr_parser_var->get_error() << ")" << endl;
                                    break;
                                case evaluator3::JIT_INLINE:
                                    if(!curr_parser_var->compile_inline())
                                        cout << "Error in " << __FILE__ << ":" << __LINE__
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
                            else cout << "[Config] Unsupported param \"" << param << "\" in section \"" << section
                                      << (subsection.empty() ? string("") : (string(".") + subsection)) << "\"" << endl;
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
                            else cout << "[Config] Unsupported param \"" << param << "\" in section \"" << section
                                      << (subsection.empty() ? string("") : (string(".") + subsection)) << "\"" << endl;
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

    if(filename.empty())
    {
        cout << "[PML Config] Warning: \"filename_pml\" parameter is not set" << endl;
        return true;
    }

    ifstream ifs(filename.c_str());
    if(!ifs.good())
    {
        cout << "[PML Config] Error in " << __FILE__ << ":" << __LINE__
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
                if(subsection.empty())
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
                                else cout << "[PML Config] Unsupported param \"" << param << "\" in section \"" << section
                                          << (subsection.empty() ? string("") : (string(".") + subsection)) << "\"" << endl;
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
                                else cout << "[PML Config] Unsupported param \"" << param << "\" in section \"" << section
                                          << (subsection.empty() ? string("") : (string(".") + subsection)) << "\"" << endl;
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
        pml_config_parameter * curr_par = &(it->second);
        if(curr_par->chi.real() >= big_num) curr_par->chi.real(def.chi.real());
        if(curr_par->chi.imag() >= big_num) curr_par->chi.imag(def.chi.imag());
        if(curr_par->m >= big_num)          curr_par->m = def.m;
        if(curr_par->width >= big_num)      curr_par->width = def.width;
    }

    return true;
}

// ============================================================================

tet_integration_config::tet_integration_config()
{
    init(0);
}

void tet_integration_config::set(size_t num, const double weights[], const double points[][4])
{
    gauss_num = num;
    gauss_weights.resize(gauss_num);
    gauss_points_master.resize(gauss_num, 4);
    for(size_t i = 0; i < gauss_num; i++)
    {
        gauss_weights[i] = weights[i];
        for(size_t j = 0; j < 4; j++)
            gauss_points_master[i][j] = points[i][j];
    }
}

void tet_integration_config::init(size_t order)
{
    switch(order)
    {
    case 2:
    {
        namespace curr = tet_integration_2::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 3:
    {
        namespace curr = tet_integration_3::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 4:
    {
        namespace curr = tet_integration_4::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 5:
    {
        namespace curr = tet_integration_5::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 6:
    {
        namespace curr = tet_integration_6::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 7:
    {
        namespace curr = tet_integration_7::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 8:
    {
        namespace curr = tet_integration_8::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 9:
    {
        namespace curr = tet_integration_9::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 10:
    {
        namespace curr = tet_integration_10::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 11:
    {
        namespace curr = tet_integration_11::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 12:
    {
        namespace curr = tet_integration_12::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 13:
    {
        namespace curr = tet_integration_13::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 14:
    {
        namespace curr = tet_integration_14::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    default:
    {
        namespace curr = tet_integration_8::tet_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    }
}

// ============================================================================

tr_integration_config::tr_integration_config()
{
    init(0);
}

void tr_integration_config::set(size_t num, const double weights[], const double points[][3])
{
    gauss_num = num;
    gauss_weights.resize(gauss_num);
    gauss_points_master.resize(gauss_num, 3);
    for(size_t i = 0; i < gauss_num; i++)
    {
        gauss_weights[i] = weights[i];
        for(size_t j = 0; j < 3; j++)
            gauss_points_master[i][j] = points[i][j];
    }
}

void tr_integration_config::init(size_t order)
{
    switch(order)
    {
    case 2:
    {
        namespace curr = tr_integration_2::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 3:
    {
        namespace curr = tr_integration_3::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 4:
    case 5:
    {
        namespace curr = tr_integration_5::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 6:
    {
        namespace curr = tr_integration_6::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 7:
    {
        namespace curr = tr_integration_7::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 8:
    {
        namespace curr = tr_integration_8::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 9:
    {
        namespace curr = tr_integration_9::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 10:
    {
        namespace curr = tr_integration_10::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 11:
    {
        namespace curr = tr_integration_11::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 12:
    {
        namespace curr = tr_integration_12::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 13:
    {
        namespace curr = tr_integration_13::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 14:
    {
        namespace curr = tr_integration_14::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 15:
    {
        namespace curr = tr_integration_15::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 16:
    {
        namespace curr = tr_integration_16::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 17:
    {
        namespace curr = tr_integration_17::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 18:
    {
        namespace curr = tr_integration_18::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 19:
    {
        namespace curr = tr_integration_19::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 20:
    {
        namespace curr = tr_integration_20::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 21:
    {
        namespace curr = tr_integration_21::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 22:
    {
        namespace curr = tr_integration_22::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 23:
    {
        namespace curr = tr_integration_23::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 24:
    {
        namespace curr = tr_integration_24::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 25:
    {
        namespace curr = tr_integration_25::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 26:
    {
        namespace curr = tr_integration_26::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 27:
    {
        namespace curr = tr_integration_27::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 28:
    {
        namespace curr = tr_integration_28::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    case 29:
    {
        namespace curr = tr_integration_29::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    default:
    {
        namespace curr = tr_integration_8::tr_integration;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master);
        break;
    }
    }
}
