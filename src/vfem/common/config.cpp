#include "config.h"

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
    v_cycle_enabled = false;
    max_iter = 15000;
    max_iter_v_cycle_local = 500;
#if !defined(USE_MKL)
    solver_name = "COCG_LLT_Smooth";
    solver_name_bound = "COCG_LLT_Smooth";
    solver_name_v_cycle_full = "COCG_LLT_Smooth";
    solver_name_v_cycle_ker = "COCG_LLT_Smooth";
#else
    solver_name = "COCG_LLT_Smooth_MKL";
    solver_name_bound = "COCG_LLT_Smooth_MKL";
    solver_name_v_cycle_full = "COCG_LLT_Smooth_MKL";
    solver_name_v_cycle_ker = "COCG_LLT_Smooth_MKL";
#endif
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

    list<string> whitelist;
    whitelist.push_back("VFEM");
    whitelist.push_back("Boundary");
    whitelist.push_back("Right");
    whitelist.push_back("Analytical");
    whitelist.push_back("Postprocessing");
    if(!cfg_file.check_sections(whitelist))
    {
        cout << "[Config] Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return init(false);
    }
    whitelist.clear();

    // Секция VFEM

    whitelist.push_back("basis_order");
    whitelist.push_back("basis_type");
    whitelist.push_back("eps_slae");
    whitelist.push_back("eps_slae_bound");
    whitelist.push_back("gamma_v_cycle_full");
    whitelist.push_back("gamma_v_cycle_ker");
    whitelist.push_back("gamma_v_cycle_0");
    whitelist.push_back("max_iter_v_cycle_local");
    whitelist.push_back("max_iter");
    whitelist.push_back("solver_name");
    whitelist.push_back("solver_name_bound");
    whitelist.push_back("solver_name_v_cycle_full");
    whitelist.push_back("solver_name_v_cycle_ker");
    whitelist.push_back("filename_mesh");
    whitelist.push_back("filename_phys");
    whitelist.push_back("filename_slae");
    whitelist.push_back("filename_pml");
    whitelist.push_back("jit_type");
    whitelist.push_back("v_cycle_enabled");
    whitelist.push_back("tet_integration_order");
    whitelist.push_back("tr_integration_order");
    if(!cfg_file.check_parameters("VFEM", whitelist))
    {
        cout << "[Config] Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return init(false);
    }
    whitelist.clear();

    map<evaluator3::jit_types, string> jit_types_table;
    jit_types_table[evaluator3::JIT_DISABLE] = "disable";
    jit_types_table[evaluator3::JIT_INLINE]  = "inline";
    jit_types_table[evaluator3::JIT_EXTCALL] = "extcall";

    basis.order              = cfg_file.get("vfem", "", "basis_order",              basis.order);
    basis.type               = cfg_file.get("vfem", "", "basis_type",               basis.type);
    eps_slae                 = cfg_file.get("vfem", "", "eps_slae",                 eps_slae);
    eps_slae_bound           = cfg_file.get("vfem", "", "eps_slae_bound",           eps_slae_bound);
    gamma_v_cycle_full       = cfg_file.get("vfem", "", "gamma_v_cycle_full",       gamma_v_cycle_full);
    gamma_v_cycle_ker        = cfg_file.get("vfem", "", "gamma_v_cycle_ker",        gamma_v_cycle_ker);
    gamma_v_cycle_0          = cfg_file.get("vfem", "", "gamma_v_cycle_0",          gamma_v_cycle_0);
    max_iter_v_cycle_local   = cfg_file.get("vfem", "", "max_iter_v_cycle_local",   max_iter_v_cycle_local);
    max_iter                 = cfg_file.get("vfem", "", "max_iter",                 max_iter);
    solver_name              = cfg_file.get("vfem", "", "solver_name",              solver_name);
    solver_name_bound        = cfg_file.get("vfem", "", "solver_name_bound",        solver_name_bound);
    solver_name_v_cycle_full = cfg_file.get("vfem", "", "solver_name_v_cycle_full", solver_name_v_cycle_full);
    solver_name_v_cycle_ker  = cfg_file.get("vfem", "", "solver_name_v_cycle_ker",  solver_name_v_cycle_ker);
    filename_mesh            = cfg_file.get("vfem", "", "filename_mesh",            filename_mesh);
    filename_phys            = cfg_file.get("vfem", "", "filename_phys",            filename_phys);
    filename_slae            = cfg_file.get("vfem", "", "filename_slae",            filename_slae);
    filename_pml             = cfg_file.get("vfem", "", "filename_pml",             filename_pml);
    string jit_type_str      = cfg_file.get("vfem", "", "jit_type",                 jit_types_table[jit_type]);
    v_cycle_enabled          = cfg_file.get("vfem", "", "v_cycle_enabled",          v_cycle_enabled);

    int tet_order            = cfg_file.get("vfem", "", "tet_integration_order",    (int)-1);
    if(tet_order > 0)
        basis.tet_int.init((size_t)tet_order);

    int tr_order             = cfg_file.get("vfem", "", "tr_integration_order",     (int)-1);
    if(tr_order > 0)
        basis.tr_int.init((size_t)tr_order);

    jit_type_str = to_lowercase(jit_type_str);
    for(map<evaluator3::jit_types, string>::iterator it = jit_types_table.begin(); it != jit_types_table.end(); ++it)
    {
        if(it->second == jit_type_str)
        {
            jit_type = it->first;
            break;
        }
    }

    // Секции Boundary, Right и Analytical

    whitelist.push_back("x");
    whitelist.push_back("y");
    whitelist.push_back("z");
    whitelist.push_back("enabled");
    if(!cfg_file.check_parameters("Boundary", whitelist) ||
       !cfg_file.check_parameters("Right", whitelist) ||
       !cfg_file.check_parameters("Analytical", whitelist))
    {
        cout << "[Config] Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return init(false);
    }
    whitelist.clear();

    trio<string, evaluator3 *, bool *> evaluators_and_sections[3] =
    {
        make_trio(string("Boundary"),   & boundary,     & boundary_enabled),
        make_trio(string("Right"),      & right,        & right_enabled),
        make_trio(string("Analytical"), & analytical,   & analytical_enabled)
    };

    for(size_t k = 0; k < 3; k++)
    {
        *(evaluators_and_sections[k].third) = cfg_file.get(evaluators_and_sections[k].first,
                "", "enabled", *(evaluators_and_sections[k].third));
        // Если секция отключена, то рассматривать дальнейшие параметры нет смысла
        if(*(evaluators_and_sections[k].third))
        {
            // Умолчательные значения следует взять из пустой подсекции
            string def_values[3] =
            {
                cfg_file.get(evaluators_and_sections[k].first, "", "x", string("0.0")),
                cfg_file.get(evaluators_and_sections[k].first, "", "y", string("0.0")),
                cfg_file.get(evaluators_and_sections[k].first, "", "z", string("0.0"))
            };
            // А теперь уже можно пройтись по всем подсекциям
            list<size_t> subsections = cfg_file.enumerate(evaluators_and_sections[k].first, (size_t*)NULL);
            for(list<size_t>::iterator it = subsections.begin(); it != subsections.end(); ++it)
            {
                size_t subsection = * it;
                array_t<evaluator_helmholtz, 3> * curr_parser = NULL;
                if(subsection == 0)
                    curr_parser = & evaluators_and_sections[k].second->default_value;
                else
                    curr_parser = & evaluators_and_sections[k].second->values[subsection];
                string values[3] =
                {
                    cfg_file.get(evaluators_and_sections[k].first, subsection, "x", string(def_values[0])),
                    cfg_file.get(evaluators_and_sections[k].first, subsection, "y", string(def_values[1])),
                    cfg_file.get(evaluators_and_sections[k].first, subsection, "z", string(def_values[2]))
                };
                for(size_t i = 0; i < 3; i++)
                {
                    if(!((*curr_parser)[i].parse(values[i]) && (*curr_parser)[i].simplify()))
                    {
                        cout << "Error in " << __FILE__ << ":" << __LINE__
                             << " while parsing " << values[i] << " ("
                             << (*curr_parser)[i].get_error() << ")" << endl;
                        return init(false);
                    }
                    switch(jit_type)
                    {
                    case evaluator3::JIT_EXTCALL:
                        if(!(*curr_parser)[i].compile_extcall())
                            cout << "Error in " << __FILE__ << ":" << __LINE__
                                 << " while compiling (extcall) " << values[i] << " ("
                                 << (*curr_parser)[i].get_error() << ")" << endl;
                        break;
                    case evaluator3::JIT_INLINE:
                        if(!(*curr_parser)[i].compile_inline())
                            cout << "Error in " << __FILE__ << ":" << __LINE__
                                 << " while compiling (inline) " << values[i] << " ("
                                 << (*curr_parser)[i].get_error() << ")" << endl;
                        break;
                    default:
                        break;
                    }
                }
            }
        }
    }

    // Секция Postprocessing

    whitelist.push_back("filename");
    whitelist.push_back("slice_var_name");
    whitelist.push_back("slice_var_value");
    whitelist.push_back("var1_name");
    whitelist.push_back("var1_from");
    whitelist.push_back("var1_to");
    whitelist.push_back("var1_num");
    whitelist.push_back("var2_name");
    whitelist.push_back("var2_from");
    whitelist.push_back("var2_to");
    whitelist.push_back("var2_num");
    whitelist.push_back("line_var1_name");
    whitelist.push_back("line_var1_value");
    whitelist.push_back("line_var2_name");
    whitelist.push_back("line_var2_value");
    whitelist.push_back("var_name");
    whitelist.push_back("var_from");
    whitelist.push_back("var_to");
    whitelist.push_back("var_num");
    whitelist.push_back("timestamp");
    whitelist.push_back("type");
    if(!cfg_file.check_parameters("Postprocessing", whitelist))
    {
        cout << "[Config] Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return init(false);
    }
    whitelist.clear();

    map<unsigned char, string> post_types_table;
    post_types_table[1] = "1d";
    post_types_table[2] = "2d";
    post_types_table[3] = "3d";

    list<size_t> post_subsections = cfg_file.enumerate("Postprocessing", (size_t*)NULL);
    for(list<size_t>::iterator it = post_subsections.begin(); it != post_subsections.end(); ++it)
    {
        size_t subsection = * it;
        postprocessor * p = &(post[subsection]);
        p->filename                 = cfg_file.get("postprocessing", subsection, "filename",        p->filename);
        p->param_2d.slice_var_name  = cfg_file.get("postprocessing", subsection, "slice_var_name",  p->param_2d.slice_var_name);
        p->param_2d.slice_var_value = cfg_file.get("postprocessing", subsection, "slice_var_value", p->param_2d.slice_var_value);
        p->param_2d.var1_name       = cfg_file.get("postprocessing", subsection, "var1_name",       p->param_2d.var1_name);
        p->param_2d.var1_from       = cfg_file.get("postprocessing", subsection, "var1_from",       p->param_2d.var1_from);
        p->param_2d.var1_to         = cfg_file.get("postprocessing", subsection, "var1_to",         p->param_2d.var1_to);
        p->param_2d.var1_num        = cfg_file.get("postprocessing", subsection, "var1_num",        p->param_2d.var1_num);
        p->param_2d.var2_name       = cfg_file.get("postprocessing", subsection, "var2_name",       p->param_2d.var2_name);
        p->param_2d.var2_from       = cfg_file.get("postprocessing", subsection, "var2_from",       p->param_2d.var2_from);
        p->param_2d.var2_to         = cfg_file.get("postprocessing", subsection, "var2_to",         p->param_2d.var2_to);
        p->param_2d.var2_num        = cfg_file.get("postprocessing", subsection, "var2_num",        p->param_2d.var2_num);
        p->param_1d.line_var1_name  = cfg_file.get("postprocessing", subsection, "line_var1_name",  p->param_1d.line_var1_name);
        p->param_1d.line_var1_value = cfg_file.get("postprocessing", subsection, "line_var1_value", p->param_1d.line_var1_value);
        p->param_1d.line_var2_name  = cfg_file.get("postprocessing", subsection, "line_var2_name",  p->param_1d.line_var2_name);
        p->param_1d.line_var2_value = cfg_file.get("postprocessing", subsection, "line_var2_value", p->param_1d.line_var2_value);
        p->param_1d.var_name        = cfg_file.get("postprocessing", subsection, "var_name",        p->param_1d.var_name);
        p->param_1d.var_from        = cfg_file.get("postprocessing", subsection, "var_from",        p->param_1d.var_from);
        p->param_1d.var_to          = cfg_file.get("postprocessing", subsection, "var_to",          p->param_1d.var_to);
        p->param_1d.var_num         = cfg_file.get("postprocessing", subsection, "var_num",         p->param_1d.var_num);
        p->timestamp                = cfg_file.get("postprocessing", subsection, "timestamp",       p->timestamp);
        string post_type            = cfg_file.get("postprocessing", subsection, "type",            post_types_table[p->type]);
        for(map<unsigned char, string>::const_iterator it = post_types_table.begin(); it != post_types_table.end(); ++it)
        {
            if(it->second == post_type)
            {
                p->type = it->first;
                break;
            }
        }
    }
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

    inifile cfg_file(filename);
    if(!cfg_file.good())
    {
        cout << "[PML Config] Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return false;
    }

    list<string> whitelist;
    whitelist.push_back("PML");
    if(!cfg_file.check_sections(whitelist))
    {
        cout << "[PML Config] Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return init(false);
    }
    whitelist.clear();

    // Секция PML

    whitelist.push_back("chi_imag");
    whitelist.push_back("chi_real");
    whitelist.push_back("m");
    whitelist.push_back("width");
    whitelist.push_back("x0");
    whitelist.push_back("x1");
    whitelist.push_back("y0");
    whitelist.push_back("y1");
    whitelist.push_back("z0");
    whitelist.push_back("z1");
    if(!cfg_file.check_parameters("PML", whitelist))
    {
        cout << "[PML Config] Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << filename << endl;
        return init(false);
    }
    whitelist.clear();

    // Умолчательные значения зададим пока огромными числами
    complex<double> chi_def(cfg_file.get("PML", "", "chi_real", DBL_MAX), cfg_file.get("PML", "", "chi_imag", DBL_MAX));
    double m_def     = cfg_file.get("PML", "", "m",     DBL_MAX);
    double width_def = cfg_file.get("PML", "", "width", DBL_MAX);
    phys_pml.x0      = cfg_file.get("PML", "", "x0",    DBL_MAX);
    phys_pml.x1      = cfg_file.get("PML", "", "x1",   -DBL_MAX);
    phys_pml.y0      = cfg_file.get("PML", "", "y0",    DBL_MAX);
    phys_pml.y1      = cfg_file.get("PML", "", "y1",   -DBL_MAX);
    phys_pml.z0      = cfg_file.get("PML", "", "z0",    DBL_MAX);
    phys_pml.z1      = cfg_file.get("PML", "", "z1",   -DBL_MAX);

    list<size_t> pml_subsections = cfg_file.enumerate("PML", (size_t*)NULL);
    for(list<size_t>::iterator it = pml_subsections.begin(), it_end = pml_subsections.end(); it != it_end; ++it)
    {
        size_t pml_phys = * it;
        if(pml_phys != 0)
        {
            pml_config_parameter * curr = & (phys_pml.params[pml_phys]);
            double chi_re = cfg_file.get("PML", pml_phys, "chi_real", chi_def.real());
            double chi_im = cfg_file.get("PML", pml_phys, "chi_imag", chi_def.imag());
            curr->chi = complex<double>(chi_re, chi_im);
            curr->m     = cfg_file.get("PML", pml_phys, "m",     m_def);
            curr->width = cfg_file.get("PML", pml_phys, "width", width_def);
        }
    }

    // Проверим введенное на адекватность
    pml_config_parameter def;
    double big_num = DBL_MAX * 0.999;
    if(chi_def.real() < big_num) def.chi = complex<double>(chi_def.real(), def.chi.imag());
    if(chi_def.imag() < big_num) def.chi = complex<double>(def.chi.real(), chi_def.imag());
    if(m_def < big_num)          def.m = m_def;
    if(width_def < big_num)      def.width = width_def;
    for(map<size_t, pml_config_parameter>::iterator it = phys_pml.params.begin(); it != phys_pml.params.end(); ++it)
    {
        pml_config_parameter * curr_par = &(it->second);
        if(curr_par->chi.real() >= big_num) curr_par->chi = complex<double>(def.chi.real(), curr_par->chi.imag());
        if(curr_par->chi.imag() >= big_num) curr_par->chi = complex<double>(curr_par->chi.real(), def.chi.imag());
        if(curr_par->m >= big_num)          curr_par->m = def.m;
        if(curr_par->width >= big_num)      curr_par->width = def.width;
    }

    return true;
}
