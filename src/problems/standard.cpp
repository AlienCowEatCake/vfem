#include "problems.h"

#if defined PROBLEM_STANDARD
#if defined VFEM_USE_PML
#error "Please, reconfigure!"
#endif

cvector3 func_true(const config_type * config, const point & p, const phys_area & phys)
{
    cvector3 result(0, 0, 0);
    if(!config->analytical_enabled) return result;
    complex<double> k2(- phys.omega * phys.omega * phys.epsilon, phys.omega * phys.sigma);
    config_type * cfg = const_cast<config_type *>(config);
    map<size_t, array_t<parser<complex<double> >, 3> >::iterator it = cfg->analytical.values.find(phys.gmsh_num);
    parser<complex<double> > * evaluator[3];
    if(it == cfg->analytical.values.end())
        for(size_t i = 0; i < 3; i++)
        {
            evaluator[i] = &(cfg->analytical.default_value[i]);
            evaluator[i]->set_const("epsilon", phys.epsilon);
            evaluator[i]->set_const("sigma", phys.sigma);
            evaluator[i]->set_const("mu", phys.mu);
            evaluator[i]->set_const("J0", phys.J0);
            evaluator[i]->set_const("k2", k2);
        }
    else
        for(size_t i = 0; i < 3; i++)
            evaluator[i] = &(it->second[i]);
    for(size_t i = 0; i < 3; i++)
    {
        evaluator[i]->set_const("x", p.x);
        evaluator[i]->set_const("y", p.y);
        evaluator[i]->set_const("z", p.z);
        bool status = evaluator[i]->calculate(result[i]);
        if(!status) cerr << "[Parser] " << evaluator[i]->get_error() << endl;
    }
    return result;
}

cvector3 func_rp(const config_type * config, const point & p, const phys_area & phys)
{
    cvector3 result(0, 0, 0);
    if(!config->right_enabled) return result;
    complex<double> k2(- phys.omega * phys.omega * phys.epsilon, phys.omega * phys.sigma);
    config_type * cfg = const_cast<config_type *>(config);
    map<size_t, array_t<parser<complex<double> >, 3> >::iterator it = cfg->right.values.find(phys.gmsh_num);
    parser<complex<double> > * evaluator[3];
    if(it == cfg->right.values.end())
        for(size_t i = 0; i < 3; i++)
        {
            evaluator[i] = &(cfg->right.default_value[i]);
            evaluator[i]->set_const("epsilon", phys.epsilon);
            evaluator[i]->set_const("sigma", phys.sigma);
            evaluator[i]->set_const("mu", phys.mu);
            evaluator[i]->set_const("J0", phys.J0);
            evaluator[i]->set_const("k2", k2);
        }
    else
        for(size_t i = 0; i < 3; i++)
            evaluator[i] = &(it->second[i]);
    for(size_t i = 0; i < 3; i++)
    {
        evaluator[i]->set_const("x", p.x);
        evaluator[i]->set_const("y", p.y);
        evaluator[i]->set_const("z", p.z);
        bool status = evaluator[i]->calculate(result[i]);
        if(!status) cerr << "[Parser] " << evaluator[i]->get_error() << endl;
    }
    return result;
}

cvector3 func_b1(const config_type * config, const point & p, const phys_area & phys)
{
    cvector3 result(0, 0, 0);
    if(!config->boundary_enabled) return result;
    complex<double> k2(- phys.omega * phys.omega * phys.epsilon, phys.omega * phys.sigma);
    config_type * cfg = const_cast<config_type *>(config);
    map<size_t, array_t<parser<complex<double> >, 3> >::iterator it = cfg->boundary.values.find(phys.gmsh_num);
    parser<complex<double> > * evaluator[3];
    if(it == cfg->boundary.values.end())
        for(size_t i = 0; i < 3; i++)
        {
            evaluator[i] = &(cfg->boundary.default_value[i]);
            evaluator[i]->set_const("epsilon", phys.epsilon);
            evaluator[i]->set_const("sigma", phys.sigma);
            evaluator[i]->set_const("mu", phys.mu);
            evaluator[i]->set_const("J0", phys.J0);
            evaluator[i]->set_const("k2", k2);
        }
    else
        for(size_t i = 0; i < 3; i++)
            evaluator[i] = &(it->second[i]);
    for(size_t i = 0; i < 3; i++)
    {
        evaluator[i]->set_const("x", p.x);
        evaluator[i]->set_const("y", p.y);
        evaluator[i]->set_const("z", p.z);
        bool status = evaluator[i]->calculate(result[i]);
        if(!status) cerr << "[Parser] " << evaluator[i]->get_error() << endl;
    }
    return result;
}

void postprocessing(VFEM & v, char * timebuf)
{
    MAYBE_UNUSED(v);
    MAYBE_UNUSED(timebuf);

    if(v.config.filename_slae != "")
        v.slae.dump_x(v.config.filename_slae);

    if(v.config.analytical_enabled)
        v.calculate_diff();

    for(map<size_t, postprocessor>::iterator it = v.config.post.begin(); it != v.config.post.end(); ++it)
    {
        string new_name;
        postprocessor * p = &(it->second);
        if(p->timestamp)
        {
            size_t dot_pos = p->filename.find_last_of(".");
            if(dot_pos != string::npos)
                new_name = p->filename.substr(0, dot_pos) + "_" + string(timebuf) + p->filename.substr(dot_pos);
            else
                new_name = p->filename + "_" + string(timebuf);
        }
        else
            new_name = p->filename;

        switch(p->type)
        {
        case 1:
            v.output_line(new_name, p->param_1d.line_var1_name, p->param_1d.line_var1_value,
                          p->param_1d.line_var2_name, p->param_1d.line_var2_value,
                          p->param_1d.var_name, p->param_1d.var_from, p->param_1d.var_to,
                          p->param_1d.var_num);
            break;
        case 2:
            v.output_slice(new_name, p->param_2d.slice_var_name, p->param_2d.slice_var_value,
                           p->param_2d.var1_name, p->param_2d.var1_from, p->param_2d.var1_to,
                           p->param_2d.var1_num, p->param_2d.var2_name, p->param_2d.var2_from,
                           p->param_2d.var2_to, p->param_2d.var2_num);
            break;
        case 3:
            v.output(new_name);
            break;
        }
    }
}

#endif
