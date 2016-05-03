#include "problems.h"

// Функция аналитического решения
cvector3 func_true(const point & p, const phys_area & phys, void * data)
{
    (void)(phys);
    pair<const config_type *, array_t<evaluator_helmholtz *, 3> > * params =
            (pair<const config_type *, array_t<evaluator_helmholtz *, 3> > *)(data);
    cvector3 result(0, 0, 0);
    if(!params->first->analytical_enabled) return result;

    // TODO: Костыль
    phys_area * phys_editable = const_cast<phys_area *>(&phys);
    phys_editable->sigma[threads_config::analytical].set_x(p.x);
    phys_editable->sigma[threads_config::analytical].set_y(p.y);
    phys_editable->sigma[threads_config::analytical].set_z(p.z);
    double sigma;
    if(!phys_editable->sigma[threads_config::analytical].calculate(sigma))
        cout << "[Parser] " << phys_editable->sigma[threads_config::analytical].get_error() << endl;
    complex<double> k2(- phys.epsilon * phys.omega * phys.omega, phys.omega * sigma);

    for(size_t i = 0; i < 3; i++)
    {
        evaluator_helmholtz * e = params->second[i];
        e->set_x(p.x);
        e->set_y(p.y);
        e->set_z(p.z);
        e->set_sigma(sigma);
        e->set_k2(k2);
        bool status = e->calculate(result[i]);
        if(!status) cout << "[Parser] " << e->get_error() << endl;
    }
    return result;
}

// Правая часть
cvector3 func_rp(const point & p, const phys_area & phys, void * data)
{
    (void)(phys);
    pair<const config_type *, array_t<evaluator_helmholtz *, 3> > * params =
            (pair<const config_type *, array_t<evaluator_helmholtz *, 3> > *)(data);
    cvector3 result(0, 0, 0);
    if(!params->first->right_enabled) return result;

    // TODO: Костыль
    phys_area * phys_editable = const_cast<phys_area *>(&phys);
    phys_editable->sigma[threads_config::right_part].set_x(p.x);
    phys_editable->sigma[threads_config::right_part].set_y(p.y);
    phys_editable->sigma[threads_config::right_part].set_z(p.z);
    double sigma;
    if(!phys_editable->sigma[threads_config::right_part].calculate(sigma))
        cout << "[Parser] " << phys_editable->sigma[threads_config::right_part].get_error() << endl;
    complex<double> k2(- phys.epsilon * phys.omega * phys.omega, phys.omega * sigma);

    for(size_t i = 0; i < 3; i++)
    {
        evaluator_helmholtz * e = params->second[i];
        e->set_x(p.x);
        e->set_y(p.y);
        e->set_z(p.z);
        e->set_sigma(sigma);
        e->set_k2(k2);
        bool status = e->calculate(result[i]);
        if(!status) cout << "[Parser] " << e->get_error() << endl;
    }
    return result;
}

// Функция неоднородных первых краевых условий
cvector3 func_b1(const point & p, const phys_area & phys, void * data)
{
    (void)(phys);
    pair<const config_type *, array_t<evaluator_helmholtz *, 3> > * params =
            (pair<const config_type *, array_t<evaluator_helmholtz *, 3> > *)(data);
    cvector3 result(0, 0, 0);
    if(!params->first->boundary_enabled) return result;

    // TODO: Костыль
    phys_area * phys_editable = const_cast<phys_area *>(&phys);
    phys_editable->sigma[threads_config::boundary].set_x(p.x);
    phys_editable->sigma[threads_config::boundary].set_y(p.y);
    phys_editable->sigma[threads_config::boundary].set_z(p.z);
    double sigma;
    if(!phys_editable->sigma[threads_config::boundary].calculate(sigma))
        cout << "[Parser] " << phys_editable->sigma[threads_config::boundary].get_error() << endl;
    complex<double> k2(- phys.epsilon * phys.omega * phys.omega, phys.omega * sigma);

    for(size_t i = 0; i < 3; i++)
    {
        evaluator_helmholtz * e = params->second[i];
        e->set_x(p.x);
        e->set_y(p.y);
        e->set_z(p.z);
        e->set_sigma(sigma);
        e->set_k2(k2);
        bool status = e->calculate(result[i]);
        if(!status) cout << "[Parser] " << e->get_error() << endl;
    }
    return result;
}

// Постпроцессор
void postprocessing(VFEM & v, const char * timebuf)
{
    if(!v.config.filename_slae.empty())
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
