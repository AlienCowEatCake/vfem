#include "vfem.h"

// Получение степеней свободы тетраэдра в глобальной матрице
size_t VFEM::get_tet_dof(const tetrahedron_base * fe, size_t i) const
{
    assert(i < config.basis.tet_bf_num);
    if(i < 6)       return fe->edges[i]->num;
    else if(i < 12) return fe->edges[i-6]->num + edges.size();
    else if(i < 16) return fe->faces[i-12]->num + 2 * edges.size();
    else if(i < 20) return fe->faces[i-16]->num + 2 * edges.size() + faces.size();
    else if(i < 24) return fe->faces[i-20]->num + 2 * edges.size() + 2 * faces.size();
    else if(i < 30) return fe->edges[i-24]->num + 2 * edges.size() + 3 * faces.size();
    return 0;
}

// Получение степеней свободы тетраэдра в матрице ядра
size_t VFEM::get_tet_ker_dof(const tetrahedron_base * fe, size_t i) const
{
    assert(i < config.basis.tet_ker_bf_num);
    if(i < 4)   return fe->nodes[i]->num;
    if(i < 10)  return fe->edges[i-4]->num + nodes.size();
    if(i < 14)  return fe->faces[i-10]->num + edges.size() + nodes.size();
    if(i < 20)  return fe->edges[i-14]->num + faces.size() + edges.size() + nodes.size();
    return 0;
}

// Получение степеней свободы треугольника в глобальной матрице
size_t VFEM::get_tr_dof(const triangle * tr, size_t i) const
{
    assert(i < config.basis.tr_bf_num);
    if(i < 3)       return tr->edges[i]->num;
    else if(i < 6)  return tr->edges[i-3]->num + edges.size();
    else if(i < 7)  return tr->faces->num + 2 * edges.size();
    else if(i < 8)  return tr->faces->num + 2 * edges.size() + faces.size();
    else if(i < 9)  return tr->faces->num + 2 * edges.size() + 2 * faces.size();
    else if(i < 12) return tr->edges[i-9]->num + 2 * edges.size() + 3 * faces.size();
    return 0;
}

// Получение степеней свободы треугольника в матрице ядра
size_t VFEM::get_tr_ker_dof(const triangle * tr, size_t i) const
{
    assert(i < config.basis.tr_ker_bf_num);
    if(i < 3)       return tr->nodes[i]->num;
    else if(i < 6)  return tr->edges[i-3]->num + nodes.size();
    else if(i < 7)  return tr->faces->num + edges.size() + nodes.size();
    else if(i < 10) return tr->edges[i-7]->num + edges.size() + faces.size() + nodes.size();
    return 0;
}

// Получение степеней свободы треугольника в матрице по границе
size_t VFEM::get_tr_surf_dof(const triangle * tr, size_t i) const
{
    return global_to_local.find(get_tr_dof(tr, i))->second;
}

void VFEM::generate_portrait()
{
    cout << "Generating portrait ..." << endl;

    set<size_t> * portrait = new set<size_t> [dof_num];
    for(size_t k = 0; k < fes.size(); k++)
    {
        show_progress("step 1", k, fes.size());

        array_t<size_t> dof(config.basis.tet_bf_num);
        for(size_t i = 0; i < config.basis.tet_bf_num; i++)
            dof[i] = get_tet_dof(&fes[k], i);

        for(size_t i = 0; i < config.basis.tet_bf_num; i++)
        {
            size_t a = dof[i];
            for(size_t j = 0; j < i; j++)
            {
                size_t b = dof[j];
                if(b > a)
                    portrait[b].insert(a);
                else
                    portrait[a].insert(b);
            }
        }
    }

    size_t gg_size = 0;
    for(size_t i = 0; i < dof_num; i++)
        gg_size += portrait[i].size();

    slae.alloc_all(dof_num, gg_size);

    slae.ig[0] = 0;
    slae.ig[1] = 0;
    size_t tmp = 0;
    for(size_t i = 0; i < dof_num; i++)
    {
        show_progress("step 2", i, dof_num);

        for(set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
        {
            slae.jg[tmp] = *j;
            tmp++;
        }
        slae.ig[i + 1] = slae.ig[i] + portrait[i].size();

        portrait[i].clear();
    }

    delete [] portrait;
}

void VFEM::generate_ker_portrait()
{
    cout << "Generating kernel portrait ..." << endl;

    set<size_t> * portrait = new set<size_t> [ker_dof_num];
    for(size_t k = 0; k < fes.size(); k++)
    {
        show_progress("step 1", k, fes.size());

        array_t<size_t> ker_dof(config.basis.tet_ker_bf_num);
        for(size_t i = 0; i < config.basis.tet_ker_bf_num; i++)
            ker_dof[i] = get_tet_ker_dof(&fes[k], i);

        for(size_t i = 0; i < config.basis.tet_ker_bf_num; i++)
        {
            size_t a = ker_dof[i];
            for(size_t j = 0; j < i; j++)
            {
                size_t b = ker_dof[j];
                if(b > a)
                    portrait[b].insert(a);
                else
                    portrait[a].insert(b);
            }
        }
    }

    size_t gg_size = 0;
    for(size_t i = 0; i < ker_dof_num; i++)
        gg_size += portrait[i].size();

    ker_slae.alloc_all(ker_dof_num, gg_size);

    ker_slae.ig[0] = 0;
    ker_slae.ig[1] = 0;
    size_t tmp = 0;
    for(size_t i = 0; i < ker_dof_num; i++)
    {
        show_progress("step 2", i, ker_dof_num);

        for(set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
        {
            ker_slae.jg[tmp] = *j;
            tmp++;
        }
        ker_slae.ig[i + 1] = ker_slae.ig[i] + portrait[i].size();

        portrait[i].clear();
    }

    delete [] portrait;
}

void VFEM::generate_surf_portrait()
{
    cout << "Generaing surface portrait ..." << endl;

    size_t dof_surf_num = global_to_local.size();
    set<size_t> * portrait = new set<size_t> [dof_surf_num];
    for(size_t k = 0; k < trs.size(); k++)
    {
        show_progress("step 1", k, trs.size());

        array_t<size_t> dof_surf(config.basis.tr_bf_num);
        for(size_t i = 0; i < config.basis.tr_bf_num; i++)
            dof_surf[i] = get_tr_surf_dof(trs[k], i);

        if(trs[k]->phys->type_of_bounds == 1)
        {
            for(size_t i = 0; i < config.basis.tr_bf_num; i++)
            {
                size_t a = dof_surf[i];
                for(size_t j = 0; j < i; j++)
                {
                    size_t b = dof_surf[j];
                    if(b > a)
                        portrait[b].insert(a);
                    else
                        portrait[a].insert(b);
                }
            }
        }
    }

    size_t gg_size = 0;
    for(size_t i = 0; i < dof_surf_num; i++)
        gg_size += portrait[i].size();

    surf_slae.alloc_all(dof_surf_num, gg_size);

    surf_slae.ig[0] = 0;
    surf_slae.ig[1] = 0;
    size_t tmp = 0;
    for(size_t i = 0; i < dof_surf_num; i++)
    {
        show_progress("step 2", i, dof_surf_num);

        for(set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
        {
            surf_slae.jg[tmp] = *j;
            tmp++;
        }
        surf_slae.ig[i + 1] = surf_slae.ig[i] + portrait[i].size();

        portrait[i].clear();
    }

    delete [] portrait;
}

// Добавление локальных матриц от одного КЭ в глобальную
template<typename U, typename V>
void VFEM::process_fe(const U * curr_fe, const V *)
{
    // Получение физических параметров для заданного КЭ
    phys_area ph = curr_fe->get_phys_area();
    complex<double> k2(- ph.epsilon * ph.omega * ph.omega, ph.omega * ph.sigma);

    // Инициализация параметров вычислителей для правой части
    pair<const config_type *, array_t<evaluator<complex<double> > *, 3> >
            params_object(& config, array_t<evaluator<complex<double> > *, 3>());
    if(config.right_enabled)
    {
        map<size_t, array_t<evaluator<complex<double> >, 3> >::iterator
                it = config.right.values.find(ph.gmsh_num);
        if(it != config.right.values.end())
            for(size_t i = 0; i < 3; i++)
                params_object.second[i] = &(it->second[i]);
        else
            for(size_t i = 0; i < 3; i++)
            {
                evaluator<complex<double> > * ev_curr = &(config.right.default_value[i]);
                params_object.second[i] = ev_curr;
                ev_curr->set_var("epsilon", ph.epsilon);
                ev_curr->set_var("sigma", ph.sigma);
                ev_curr->set_var("mu", ph.mu);
                ev_curr->set_var("J0", ph.J0);
                ev_curr->set_var("k2", k2);
            }
    }

    // Получение степеней свободы
    array_t<size_t> dof(config.basis.tet_bf_num);
    for(size_t i = 0; i < config.basis.tet_bf_num; i++)
        dof[i] = get_tet_dof(curr_fe, i);

    // Получение локальных матриц и правой части
    V matrix_G = curr_fe->G();
    V matrix_M = curr_fe->M();
    array_t<complex<double> > array_rp = curr_fe->rp(func_rp, &params_object);

    // Основная матрица
    for(size_t i = 0; i < config.basis.tet_bf_num; i++)
    {
        complex<double> add;
        size_t i_num = dof[i];
        for(size_t j = 0; j < config.basis.tet_bf_num; j++)
        {
            size_t j_num = dof[j];
            add = matrix_G[i][j] / ph.mu + matrix_M[i][j] * k2;
            slae.add(i_num, j_num, add);
        }
        add = array_rp[i];
        slae.rp[i_num] += add;
    }

    // Матрица ядра
    if(config.v_cycle_enabled)
    {
        array_t<size_t> ker_dof(config.basis.tet_ker_bf_num);
        for(size_t i = 0; i < config.basis.tet_ker_bf_num; i++)
            ker_dof[i] = get_tet_ker_dof(curr_fe, i);

        V matrix_K = curr_fe->K();
        for(size_t i = 0; i < config.basis.tet_ker_bf_num; i++)
        {
            size_t i_num = ker_dof[i];
            for(size_t j = 0; j < config.basis.tet_ker_bf_num; j++)
                ker_slae.add(i_num, ker_dof[j], matrix_K[i][j] * k2);
        }
    }
}

void VFEM::assemble_matrix()
{
    cout << "Assembling matrix ..." << endl;

    cout << " > Assembling matrix ..." << endl;
    // Cборка основной матрицы
    for(size_t k = 0; k < fes.size(); k++)
    {
        show_progress("", k, fes.size());
#if defined VFEM_USE_PML
        if(!is_pml(fes[k].barycenter, &fes[k], &phys_pml))
            process_fe(fes[k].to_std(), (matrix_t<double> *)NULL);
        else
#endif
            process_fe(&fes[k], (l_matrix *)NULL);
    }
}

void VFEM::applying_bound()
{
    cout << " > Applying first bound ..." << endl;

    if(config.boundary_enabled)
    {

        // Решение СЛАУ по границе
        if(global_to_local.size() > 0)
        {
            for(size_t k = 0; k < trs.size(); k++)
            {
                show_progress("building matrix", k, trs.size());

                phys_area ph = trs[k]->get_phys_area();
                if(ph.type_of_bounds == 1)
                {
                    // Получение физических параметров для текущего треугольника границы
                    complex<double> k2(- ph.epsilon * ph.omega * ph.omega, ph.omega * ph.sigma);

                    // Инициализация параметров вычислителей для первого краевого
                    pair<const config_type *, array_t<evaluator<complex<double> > *, 3> >
                            params_object(& config, array_t<evaluator<complex<double> > *, 3>());
                    map<size_t, array_t<evaluator<complex<double> >, 3> >::iterator
                            it = config.boundary.values.find(ph.gmsh_num);
                    if(it != config.boundary.values.end())
                        for(size_t i = 0; i < 3; i++)
                            params_object.second[i] = &(it->second[i]);
                    else
                        for(size_t i = 0; i < 3; i++)
                        {
                            evaluator<complex<double> > * ev_curr = &(config.boundary.default_value[i]);
                            params_object.second[i] = ev_curr;
                            ev_curr->set_var("epsilon", ph.epsilon);
                            ev_curr->set_var("sigma", ph.sigma);
                            ev_curr->set_var("mu", ph.mu);
                            ev_curr->set_var("J0", ph.J0);
                            ev_curr->set_var("k2", k2);
                        }

                    // Получение степеней свободы
                    array_t<size_t> dof_surf(config.basis.tr_bf_num);
                    for(size_t i = 0; i < config.basis.tr_bf_num; i++)
                        dof_surf[i] = get_tr_surf_dof(trs[k], i);

                    // Получение локальных матриц и правой части
                    matrix_t<double> M_surf = trs[k]->M();
                    array_t<complex<double> > b_surf = trs[k]->rp(func_b1, &params_object);

                    // Добавляем полученное в матрицу с краевыми
                    for(size_t i = 0; i < config.basis.tr_bf_num; i++)
                    {
                        for(size_t j = 0; j < i; j++)
                            surf_slae.add(dof_surf[i], dof_surf[j], M_surf[i][j]);
                        surf_slae.di[dof_surf[i]] += M_surf[i][i];
                        surf_slae.rp[dof_surf[i]] += b_surf[i];
                    }
                }
            }
            surf_slae.solve(config.eps_slae_bound, config.max_iter);
        }

        // Учет первых краевых
        for(size_t k = 0; k < slae.n; k++) 	  // Пробегаем по всей матрице
        {
            show_progress("applying", k, slae.n);

            if(global_to_local.find(k) != global_to_local.end())
            {
                complex<double> val = surf_slae.x[global_to_local[k]];
                slae.rp[k] = val;
                slae.di[k] = 1.0;
                slae.x[k] = val; // Начальное приближение сразу знаем

                for(size_t i = slae.ig[k]; i < slae.ig[k + 1]; i++)
                {
//                    if(global_to_local.find(slae.jg[i]) == global_to_local.end())
//                        slae.rp[slae.jg[i]] -= slae.gg[i] * val;
//                    slae.gg[i] = 0.0;
                    if(global_to_local.find(slae.jg[i]) != global_to_local.end())
                        slae.gu[i] = 0.0;
                    slae.gl[i] = 0.0;
                }
            }
            else
            {
                for(size_t i = slae.ig[k]; i < slae.ig[k + 1]; i++)
                {
//                    if(global_to_local.find(slae.jg[i]) != global_to_local.end())
//                    {
//                        complex<double> val = surf_slae.x[global_to_local[slae.jg[i]]];
//                        slae.rp[k] -= slae.gg[i] * val;
//                        slae.gg[i] = 0.0;
//                    }
                    if(global_to_local.find(slae.jg[i]) != global_to_local.end())
                        slae.gu[i] = 0.0;
                }
            }
        }

    }
    else
    {

        for(size_t k = 0; k < slae.n; k++) 	  // Пробегаем по всей матрице
        {
            show_progress("", k, slae.n);

            if(dof_first.find(k) != dof_first.end())
            {
                slae.rp[k] = 0.0;
                slae.di[k] = 1.0;

//                for(size_t i = slae.ig[k]; i < slae.ig[k + 1]; i++)
//                    slae.gg[i] = 0.0;
                for(size_t i = slae.ig[k]; i < slae.ig[k + 1]; i++)
                {
                    if(dof_first.find(slae.jg[i]) != dof_first.end())
                        slae.gu[i] = 0.0;
                    slae.gl[i] = 0.0;
                }
            }
            else
            {
                for(size_t i = slae.ig[k]; i < slae.ig[k + 1]; i++)
                {
//                    if(dof_first.find(slae.jg[i]) != dof_first.end())
//                        slae.gg[i] = 0.0;
                    if(dof_first.find(slae.jg[i]) != dof_first.end())
                        slae.gu[i] = 0.0;
                }
            }
        }

    }

    if(config.v_cycle_enabled)
    {
        // Первые краевые для матрицы ядра
        for(size_t k = 0; k < ker_slae.n; k++)
        {
            // Пробегаем по всей матрице
            if(ker_dof_first.find(k) != ker_dof_first.end())
            {
                ker_slae.di[k] = 1.0;
//                for(size_t i = ker_slae.ig[k]; i < ker_slae.ig[k + 1]; i++)
//                    ker_slae.gg[i] = 0.0;
                for(size_t i = slae.ig[k]; i < slae.ig[k + 1]; i++)
                {
                    if(dof_first.find(slae.jg[i]) != dof_first.end())
                        ker_slae.gu[i] = 0.0;
                    ker_slae.gl[i] = 0.0;
                }
            }
            else
            {
                for(size_t i = ker_slae.ig[k]; i < ker_slae.ig[k + 1]; i++)
                {
//                    if(ker_dof_first.find(ker_slae.jg[i]) != ker_dof_first.end())
//                        ker_slae.gg[i] = 0.0;
                    if(dof_first.find(slae.jg[i]) != dof_first.end())
                        ker_slae.gu[i] = 0.0;
                }
            }
        }
    }
}

void VFEM::apply_point_sources()
{
    cout << " > Applying point sources ..." << endl;
    for(size_t k = 0; k < pss.size(); k++)
    {
        show_progress("", k, pss.size());
        finite_element * fe = get_fe(pss[k].first);
        for(size_t i = 0; i < config.basis.tet_bf_num; i++)
            slae.rp[get_tet_dof(fe, i)] += complex<double>(0.0, -1.0) * fe->phys->omega * pss[k].second * fe->w(i, pss[k].first);
    }
}

void VFEM::apply_edges_sources()
{
    cout << " > Applying edges sources ..." << endl;
    for(size_t k = 0; k < edges_src.size(); k++)
    {
        show_progress("", k, edges_src.size());
        size_t pos = edges_src[k].num; // Заносим только в роторные функции!
        slae.rp[pos] += complex<double>(0.0, -1.0) * edges_src[k].phys->J0 * edges_src[k].phys->omega * edges_src[k].direction;
    }
}

finite_element * VFEM::get_fe(const point & p) const
{
    finite_element * fe = NULL;
    fe = tree.find(p.x, p.y, p.z);
    if(!fe)
        cout << "Warning: Point " << p << " is outside area!" << endl;
    return fe;
}

cvector3 VFEM::solution(const point & p) const
{
    return solution(p, get_fe(p));
}

cvector3 VFEM::solution(const point & p, const finite_element * fe) const
{
    cvector3 result;
    if(fe)
    {
        for(size_t i = 0; i < config.basis.tet_bf_num; i++)
            result = result + slae.x[get_tet_dof(fe, i)] * fe->w(i, p);
    }
    return result;
}

cvector3 VFEM::rotor(const point & p) const
{
    return rotor(p, get_fe(p));
}

cvector3 VFEM::rotor(const point & p, const finite_element * fe) const
{
    cvector3 result;
    if(fe)
    {
        for(size_t i = 0; i < config.basis.tet_bf_num; i++)
            result = result + slae.x[get_tet_dof(fe, i)] * fe->rotw(i, p);
    }
    return result;
}

void VFEM::make_struct()
{
    if(config.boundary_enabled && global_to_local.size() > 0)
        generate_surf_portrait();
    generate_portrait();
    if(config.v_cycle_enabled)
        generate_ker_portrait();
}

void VFEM::make_data()
{
    assemble_matrix();
    apply_point_sources();
    apply_edges_sources();
    applying_bound();
}

void VFEM::calculate_diff()
{
    double diff = 0.0, norm = 0.0;
    vector3 diff_v3(0.0, 0.0, 0.0), norm_v3(0.0, 0.0, 0.0);
    cvector3 diff_cv3(0.0, 0.0, 0.0), norm_cv3(0.0, 0.0, 0.0);
    for(size_t k = 0; k < fes.size(); k++)
    {
        // Получение физических параметров для заданного КЭ
        phys_area ph = fes[k].get_phys_area();
        complex<double> k2(- ph.epsilon * ph.omega * ph.omega, ph.omega * ph.sigma);

        // Инициализация параметров вычислителей для аналитики
        pair<const config_type *, array_t<evaluator<complex<double> > *, 3> >
                params_object(& config, array_t<evaluator<complex<double> > *, 3>());
        map<size_t, array_t<evaluator<complex<double> >, 3> >::iterator
                it = config.analytical.values.find(ph.gmsh_num);
        if(it != config.analytical.values.end())
            for(size_t i = 0; i < 3; i++)
                params_object.second[i] = &(it->second[i]);
        else
            for(size_t i = 0; i < 3; i++)
            {
                evaluator<complex<double> > * ev_curr = &(config.analytical.default_value[i]);
                params_object.second[i] = ev_curr;
                ev_curr->set_var("epsilon", ph.epsilon);
                ev_curr->set_var("sigma", ph.sigma);
                ev_curr->set_var("mu", ph.mu);
                ev_curr->set_var("J0", ph.J0);
                ev_curr->set_var("k2", k2);
            }

        // Находим локальные веса, связанные с этим КЭ
        array_t<complex<double> > q_loc(config.basis.tet_bf_num);
        for(size_t i = 0; i < config.basis.tet_bf_num; i++)
            q_loc[i] = slae.x[get_tet_dof(&fes[k], i)];

        // И считаем разность в норме L2
        trio<double, vector3, cvector3> diff_tmp = fes[k].diff_normL2(q_loc, func_true, &params_object);
        trio<double, vector3, cvector3> norm_tmp = fes[k].normL2(func_true, &params_object);
        diff += diff_tmp.first;
        diff_v3 += diff_tmp.second;
        diff_cv3 += diff_tmp.third;
        norm += norm_tmp.first;
        norm_v3 += norm_tmp.second;
        norm_cv3 += norm_tmp.third;
    }
    cout << "Diff (L2):        " << sqrt(diff / norm) << endl;
    cout << "Diff (L2) [Ex]:   " << sqrt(diff_v3.x / norm_v3.x) << endl;
    cout << "Diff (L2) [Ey]:   " << sqrt(diff_v3.y / norm_v3.y) << endl;
    cout << "Diff (L2) [Ez]:   " << sqrt(diff_v3.z / norm_v3.z) << endl;
    cout << "Diff (L2) [ExR]:  " << sqrt(diff_cv3.x.real() / norm_cv3.x.real()) << endl;
    cout << "Diff (L2) [EyR]:  " << sqrt(diff_cv3.y.real() / norm_cv3.y.real()) << endl;
    cout << "Diff (L2) [EzR]:  " << sqrt(diff_cv3.z.real() / norm_cv3.z.real()) << endl;
    cout << "Diff (L2) [ExI]:  " << sqrt(diff_cv3.x.imag() / norm_cv3.x.imag()) << endl;
    cout << "Diff (L2) [EyI]:  " << sqrt(diff_cv3.y.imag() / norm_cv3.y.imag()) << endl;
    cout << "Diff (L2) [EzI]:  " << sqrt(diff_cv3.z.imag() / norm_cv3.z.imag()) << endl;
}
