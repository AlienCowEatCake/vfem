#include "vfem.h"

void VFEM::generate_portrait()
{
    cout << "Generating portrait ..." << endl;

    set<size_t> * portrait = new set<size_t> [dof_num];
    for(size_t k = 0; k < fes.size(); k++)
    {
        show_progress("step 1", k, fes.size());

        for(size_t i = 0; i < basis::tet_bf_num; i++)
        {
            size_t a = fes[k].dof[i];
            for(size_t j = 0; j < i; j++)
            {
                size_t b = fes[k].dof[j];
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

// Генерация портрета ядра для V-цикла
void VFEM::generate_ker_portrait()
{
    cout << "Generating kernel portrait ..." << endl;

    size_t n_size = nodes.size() + edges.size();

    set<size_t> * portrait = new set<size_t> [n_size];
    for(size_t k = 0; k < fes.size(); k++)
    {
        show_progress("step 1", k, fes.size());

        // Заполняем степени свободы
        size_t dof[10];
        for(size_t i = 0; i < 4; i++)
            dof[i] = fes[k].get_node(i).num;
        for(size_t i = 0; i < 6; i++)
            dof[i + 4] = fes[k].get_edge(i).num + nodes.size();

        for(size_t i = 0; i < 10; i++)
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
    for(size_t i = 0; i < n_size; i++)
        gg_size += portrait[i].size();

    slae.ker_alloc_all(n_size, gg_size);

    slae.ker_ig[0] = 0;
    slae.ker_ig[1] = 0;
    size_t tmp = 0;
    for(size_t i = 0; i < n_size; i++)
    {
        show_progress("step 2", i, n_size);

        for(set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); j++)
        {
            slae.ker_jg[tmp] = *j;
            tmp++;
        }
        slae.ker_ig[i + 1] = slae.ker_ig[i] + portrait[i].size();

        portrait[i].clear();
    }

    delete [] portrait;
}

// Генерация портрета проектора для V-цикла
void VFEM::generate_proj_portrait()
{
    // для R
    cout << "Generating projector portrait ..." << endl;

    size_t n_size = slae.ker_n;

    set<size_t> * portrait = new set<size_t> [n_size];
    for(size_t k = 0; k < fes.size(); k++)
    {
        show_progress("step 1", k, fes.size());

        // Заполняем степени свободы
        size_t dof_lvl_2[10], dof_lvl_1[12];
        for(size_t i = 0; i < 4; i++)
            dof_lvl_2[i] = fes[k].get_node(i).num;
        for(size_t i = 0; i < 6; i++)
        {
            dof_lvl_2[i + 4] = fes[k].get_edge(i).num + nodes.size();
            dof_lvl_1[i] = fes[k].get_edge(i).num;
            dof_lvl_1[i + 6] = fes[k].get_edge(i).num + edges.size();
        }

        for(size_t i = 0; i < 10; i++)
            for(size_t j = 0; j < 12; j++)
                portrait[dof_lvl_2[i]].insert(dof_lvl_1[j]);
    }

    size_t gg_size = 0;
    for(size_t i = 0; i < n_size; i++)
        gg_size += portrait[i].size();

    slae.proj_alloc_all(n_size, gg_size);

    slae.proj_ig[0] = 0;
    slae.proj_ig[1] = 0;
    size_t tmp = 0;
    for(size_t i = 0; i < n_size; i++)
    {
        show_progress("step 2", i, n_size);

        for(set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); j++)
        {
            slae.proj_jg[tmp] = *j;
            tmp++;
        }
        slae.proj_ig[i + 1] = slae.proj_ig[i] + portrait[i].size();

        portrait[i].clear();
    }

    delete [] portrait;
}

#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
void VFEM::generate_surf_portrait()
{
    cout << "Generaing surface portrait ..." << endl;

    size_t dof_surf_num = global_to_local.size();
    set<size_t> * portrait = new set<size_t> [dof_surf_num];
    for(size_t k = 0; k < trs.size(); k++)
    {
        show_progress("step 1", k, trs.size());

        if(trs[k].phys->type_of_bounds == 1)
        {
            for(size_t i = 0; i < basis::tr_bf_num; i++)
            {
                size_t a = trs[k].dof_surf[i];
                for(size_t j = 0; j < i; j++)
                {
                    size_t b = trs[k].dof_surf[j];
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
#endif

void VFEM::assemble_matrix()
{
    cout << "Assembling matrix ..." << endl;

    cout << " > Assembling matrix ..." << endl;
    // Cборка основной матрицы
    for(size_t k = 0; k < fes.size(); k++)
    {
        show_progress("", k, fes.size());

        l_matrix matrix_G = fes[k].G();
        l_matrix matrix_M = fes[k].M();
        phys_area ph = fes[k].get_phys_area();
        array_t<complex<double>, basis::tet_bf_num> array_rp = fes[k].rp(func_rp);
        complex<double> k2(- ph.epsilon * ph.omega * ph.omega, ph.omega * ph.sigma);

        for(size_t i = 0; i < basis::tet_bf_num; i++)
        {
            complex<double> add;
            size_t i_num = fes[k].dof[i];
            for(size_t j = 0; j < i; j++)
            {
                size_t j_num = fes[k].dof[j];
                add = matrix_G[i][j] / ph.mu + matrix_M[i][j] * k2;
                slae.add(i_num, j_num, add);
            }
            add = matrix_G[i][i] / ph.mu + matrix_M[i][i] * k2;
            slae.di[i_num] += add;
            add = array_rp[i];
            slae.rp[i_num] += add;
        }

        // Степени свободы ядра
        size_t ker_dof[10];
        for(size_t i = 0; i < 4; i++)
            ker_dof[i] = fes[k].get_node(i).num;
        for(size_t i = 0; i < 6; i++)
            ker_dof[i + 4] = fes[k].get_edge(i).num + nodes.size();

        // Матрица ядра
        matrix_t<double, 10, 10> matrix_K = fes[k].K();
        for(size_t i = 0; i < 10; i++)
        {
            for(size_t j = 0; j < i; j++)
                slae.ker_add(ker_dof[i], ker_dof[j], matrix_K[i][j] * k2);
            slae.ker_di[ker_dof[i]] += matrix_K[i][i] * k2;
        }

        // Матрица проектора
        matrix_t<double, 10, 12> matrix_P = /*fes[k].GetNodalStiff();*/fes[k].P();
        for(size_t i = 0; i < 10; i++)
            for(int j = 0; j < 12; j++)
                slae.proj_set(ker_dof[i], fes[k].dof[j], matrix_P[i][j]);
    }
}

void VFEM::applying_bound()
{
    cout << " > Applying first bound ..." << endl;
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    // Решение СЛАУ по границе
    if(global_to_local.size() > 0)
    {
        for(size_t k = 0; k < trs.size(); k++)
        {
            show_progress("building matrix", k, trs.size());

            if(trs[k].phys->type_of_bounds == 1)
            {
                matrix_t<double, basis::tr_bf_num, basis::tr_bf_num> M_surf = trs[k].M();
                array_t<complex<double>, basis::tr_bf_num> b_surf = trs[k].rp(func_b1);

                for(size_t i = 0; i < basis::tr_bf_num; i++)
                {
                    for(size_t j = 0; j < i; j++)
                        surf_slae.add(trs[k].dof_surf[i], trs[k].dof_surf[j], M_surf[i][j]);
                    surf_slae.di[trs[k].dof_surf[i]] += M_surf[i][i];
                    surf_slae.rp[trs[k].dof_surf[i]] += b_surf[i];
                }
            }
        }
        extern double SLAE_SURF_EPSILON;
        surf_slae.solve(SLAE_SURF_EPSILON);
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
                if(global_to_local.find(slae.jg[i]) == global_to_local.end())
                {
                    slae.rp[slae.jg[i]] -= slae.gg[i] * val;
                }
                slae.gg[i] = 0.0;
            }
        }
        else
        {
            for(size_t i = slae.ig[k]; i < slae.ig[k + 1]; i++)
            {
                if(global_to_local.find(slae.jg[i]) != global_to_local.end())
                {
                    complex<double> val = surf_slae.x[global_to_local[slae.jg[i]]];
                    slae.rp[k] -= slae.gg[i] * val;
                    slae.gg[i] = 0.0;
                }
            }
        }
    }
#else
    for(size_t k = 0; k < slae.n; k++) 	  // Пробегаем по всей матрице
    {
        show_progress("", k, slae.n);

        if(dof_first.find(k) != dof_first.end())
        {
            slae.rp[k] = 0.0;
            slae.di[k] = 1.0;

            for(size_t i = slae.ig[k]; i < slae.ig[k + 1]; i++)
                slae.gg[i] = 0.0;
        }
        else
        {
            for(size_t i = slae.ig[k]; i < slae.ig[k + 1]; i++)
            {
                if(dof_first.find(slae.jg[i]) != dof_first.end())
                    slae.gg[i] = 0.0;
            }
        }
    }
#endif

    // Первые краевые для матрицы ядра
    for(size_t k = 0; k < slae.ker_n; k++)
    {
        // Пробегаем по всей матрице
        if(ker_edges_first.find(k) != ker_edges_first.end())
        {
            slae.ker_di[k] = 1.0;
            for(size_t i = slae.ker_ig[k]; i < slae.ker_ig[k + 1]; i++)
                slae.ker_gg[i] = 0.0;
        }
        else
        {
            for(size_t i = slae.ker_ig[k]; i < slae.ker_ig[k + 1]; i++)
            {
                if(ker_edges_first.find(slae.ker_jg[i]) != ker_edges_first.end())
                    slae.ker_gg[i] = 0.0;
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
        for(size_t i = 0; i < basis::tet_bf_num; i++)
            slae.rp[fe->dof[i]] += complex<double>(0.0, -1.0) * fe->phys->omega * pss[k].second * fe->w(i, pss[k].first);
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
        cerr << "Warning: Point " << p << " is outside area!" << endl;
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
        for(size_t i = 0; i < basis::tet_bf_num; i++)
            result = result + slae.x[fe->dof[i]] * fe->w(i, p);
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
        for(size_t i = 0; i < basis::tet_bf_num; i++)
            result = result + slae.x[fe->dof[i]] * fe->rotw(i, p);
    }
    return result;
}

void VFEM::make()
{
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    if(global_to_local.size() > 0)
        generate_surf_portrait();
#endif
    generate_portrait();
    generate_ker_portrait();
    generate_proj_portrait();
    assemble_matrix();
    apply_point_sources();
    apply_edges_sources();
    applying_bound();
}

void VFEM::solve()
{
    extern double SLAE_MAIN_EPSILON;
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    set<size_t> edges_first;
    for(map<size_t, size_t>::iterator it = global_to_local.begin(); it != global_to_local.end(); ++it)
        edges_first.insert(it->first);
    slae.solve(SLAE_MAIN_EPSILON, & edges_first, & ker_edges_first);
#else
    slae.solve(SLAE_MAIN_EPSILON, & dof_first, & ker_edges_first);
#endif
}

#if defined VFEM_USE_ANALYTICAL
void VFEM::calculate_diff() const
{
    double diff = 0.0;
    for(size_t k = 0; k < fes.size(); k++)
    {
        array_t<complex<double>, basis::tet_bf_num> q_loc;
        for(size_t i = 0; i < basis::tet_bf_num; i++)
            q_loc[i] = slae.x[fes[k].dof[i]];

        diff += fes[k].diff_normL2(q_loc, func_true);
    }
    cout << "Diff (L2): \t" << sqrt(diff) << endl;
}
#endif
