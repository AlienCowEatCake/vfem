#include "vfem.h"

VFEM::VFEM()
{
    nodes_num = 0;
    nodes = NULL;
    edges_num = 0;
    edges = NULL;
    fes_num = 0;
    fes = NULL;
    trs_num = 0;
    trs = NULL;
    bound1_num = 0;
    bound2_num = 0;
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    edges_surf_num = 0;
#endif
    pss = NULL;
    pss_num = 0;
    edges_src_num = 0;
    edges_src = NULL;
#if defined VFEM_USE_PML
    nodes_pml = NULL;
#endif
}

VFEM::~VFEM()
{
    delete [] nodes;
    delete [] edges;
    delete [] fes;
    delete [] trs;
    delete [] pss;
    delete [] edges_src;
#if defined VFEM_USE_PML
    delete [] nodes_pml;
#endif
}

void VFEM::generate_portrait()
{
    cout << "Generating portrait ..." << endl;

    size_t n_size = 2 * edges_num;

    set<size_t> * portrait = new set<size_t> [n_size];
    for(size_t k = 0; k < fes_num; k++)
    {
        show_progress("step 1", k, fes_num);

        for(size_t i = 0; i < 6; i++)
        {
            size_t a = fes[k].get_edge(i).num;
            for(size_t j = 0; j < i; j++)
            {
                size_t b = fes[k].get_edge(j).num;
                if(b > a)
                    portrait[b].insert(a);
                else
                    portrait[a].insert(b);
            }
        }
    }

    size_t gg_size = 0;

    for(size_t i = 0; i < edges_num; i++)
    {
        show_progress("step 2", i, edges_num);

        size_t ien = i + edges_num;
        portrait[ien].insert(i);
        for(set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
        {
            size_t jen = *j + edges_num;
            portrait[ien].insert(*j);
            portrait[jen].insert(i);
            portrait[ien].insert(jen);
        }
        gg_size += portrait[i].size();
    }

    for(size_t i = edges_num; i < n_size; i++)
        gg_size += portrait[i].size();

    slae.alloc_all(n_size, gg_size);

    slae.ig[0] = 0;
    slae.ig[1] = 0;
    size_t tmp = 0;
    for(size_t i = 0; i < n_size; i++)
    {
        show_progress("step 3", i, n_size);

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

#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
void VFEM::generate_surf_portrait()
{
    cout << "Generaing surface portrait ..." << endl;

    size_t n_size = 2 * edges_surf_num;

    set<size_t> * portrait = new set<size_t> [n_size];
    for(size_t k = 0; k < trs_num; k++)
    {
        show_progress("step 1", k, trs_num);

        if(trs[k].get_phys_area().type_of_bounds == 1)
        {
            for(size_t i = 0; i < 3; i++)
            {
                size_t a = trs[k].edges_surf[i];
                for(size_t j = 0; j < i; j++)
                {
                    size_t b = trs[k].edges_surf[j];
                    if(b > a)
                        portrait[b].insert(a);
                    else
                        portrait[a].insert(b);
                }
            }
        }
    }

    size_t gg_size = 0;

    for(size_t i = 0; i < edges_surf_num; i++)
    {
        show_progress("step 2", i, edges_surf_num);

        size_t ien = i + edges_surf_num;
        portrait[ien].insert(i);
        for(set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
        {
            size_t jen = *j + edges_surf_num;
            portrait[ien].insert(*j);
            portrait[jen].insert(i);
            portrait[ien].insert(jen);
        }
        gg_size += portrait[i].size();
    }

    for(size_t i = edges_surf_num; i < n_size; i++)
        gg_size += portrait[i].size();

    surf_slae.alloc_all(n_size, gg_size);

    surf_slae.ig[0] = 0;
    surf_slae.ig[1] = 0;
    size_t tmp = 0;
    for(size_t i = 0; i < n_size; i++)
    {
        show_progress("step 3", i, n_size);

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
    for(size_t k = 0; k < fes_num; k++)
    {
        show_progress("", k, fes_num);

        l_matrix12 matrix_G = fes[k].G();
        l_matrix12 matrix_M = fes[k].M();
        phys_area ph = fes[k].get_phys_area();
        carray12 array_rp = fes[k].rp(func_rp);
        complex<double> k2(- ph.epsilon * ph.omega * ph.omega, ph.omega * ph.sigma);

        for(size_t i = 0; i < 12; i++)
        {
            complex<double> add;
            size_t i_num;
            if(i < 6)   i_num = fes[k].get_edge(i).num;
            else        i_num = fes[k].get_edge(i - 6).num + edges_num;
            for(size_t j = 0; j < i; j++)
            {
                size_t j_num;
                if(j < 6)   j_num = fes[k].get_edge(j).num;
                else        j_num = fes[k].get_edge(j - 6).num + edges_num;
                add = matrix_G[i][j] / ph.mu + matrix_M[i][j] * k2;
                slae.add(i_num, j_num, add);
            }
            add = matrix_G[i][i] / ph.mu + matrix_M[i][i] * k2;
            slae.di[i_num] += add;
            add = array_rp[i];
            slae.rp[i_num] += add;
        }
    }
}

void VFEM::applying_bound()
{
    cout << " > Applying first bound ..." << endl;
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    // Решение СЛАУ по границе
    if(bound1_num > 0)
    {
        for(size_t k = 0; k < trs_num; k++)
        {
            show_progress("building matrix", k, trs_num);

            if(trs[k].get_phys_area().type_of_bounds == 1)
            {
                matrix6 M_surf = trs[k].M();
                carray6 b_surf = trs[k].rp(func_b1);
                size_t pos[6];

                for(size_t i = 0; i < 3; i++)
                {
                    pos[i] = trs[k].edges_surf[i];
                    pos[i + 3] = trs[k].edges_surf[i] + edges_surf_num;
                }

                for(size_t i = 0; i < 6; i++)
                {
                    for(size_t j = 0; j < i; j++)
                        surf_slae.add(pos[i], pos[j], M_surf[i][j]);
                    surf_slae.di[pos[i]] += M_surf[i][i];
                    surf_slae.rp[pos[i]] += b_surf[i];
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

        if(edges_first.find(k) != edges_first.end())
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
                if(edges_first.find(slae.jg[i]) != edges_first.end())
                    slae.gg[i] = 0.0;
            }
        }
    }
#endif
}

void VFEM::apply_point_sources()
{
    // Учет точечного источника
    cout << " > Applying point sources ..." << endl;
    for(size_t k = 0; k < pss_num; k++)
    {
        show_progress("", k, pss_num);
        finite_element * fe = get_fe(pss[k].first);
        for(size_t i = 0; i < 12; i++)
        {
            size_t pos;
            if(i < 6)   pos = fe->get_edge(i).num;
            else        pos = fe->get_edge(i - 6).num + edges_num;
            slae.rp[pos] += complex<double>(0.0, -1.0) * fe->get_phys_area().omega * pss[k].second * fe->w(i, pss[k].first);
        }
    }
}

void VFEM::apply_edges_sources()
{
    cout << " > Applying edges sources ..." << endl;
    for(size_t k = 0; k < edges_src_num; k++)
    {
        show_progress("", k, edges_src_num);
        size_t pos = edges_src[k].num;
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
        size_t pos[12];
        for(size_t i = 0; i < 6; i++)
        {
            pos[i] = fe->get_edge(i).num;
            pos[i + 6] = fe->get_edge(i).num + edges_num;
        }
        for(size_t i = 0; i < 12; i++)
        {
            cvector3 bfunc(fe->w(i, p));
            result = result + slae.x[pos[i]] * bfunc;
        }
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
        size_t pos[12];
        for(size_t i = 0; i < 6; i++)
        {
            pos[i] = fe->get_edge(i).num;
            pos[i + 6] = fe->get_edge(i).num + edges_num;
        }
        for(size_t i = 0; i < 12; i++)
        {
            cvector3 bfunc(fe->rotw(i, p));
            result = result + slae.x[pos[i]] * bfunc;
        }
    }
    return result;
}

void VFEM::solve()
{
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    if(bound1_num > 0)
        generate_surf_portrait();
#endif
    generate_portrait();
    assemble_matrix();
    apply_point_sources();
    apply_edges_sources();
    applying_bound();
    extern double SLAE_MAIN_EPSILON;
    slae.solve(SLAE_MAIN_EPSILON);
}

#if defined VFEM_USE_ANALYTICAL
void VFEM::calculate_diff() const
{
    double diff = 0.0;
    for(size_t k = 0; k < fes_num; k++)
    {
        size_t pos[12];
        for(size_t i = 0; i < 6; i++)
        {
            pos[i] = fes[k].get_edge(i).num;
            pos[i+6] = fes[k].get_edge(i).num + edges_num;
        }
        carray12 q_loc;
        for(int i = 0; i < 12; i++)
            q_loc[i] = slae.x[pos[i]];

        diff += fes[k].diff_normL2(q_loc, func_true);
    }
    cout << "Diff (L2): \t" << sqrt(diff) << endl;
}
#endif
