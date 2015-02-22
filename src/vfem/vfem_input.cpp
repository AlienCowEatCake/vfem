#include "vfem.h"

size_t VFEM::add_edge(edge ed, set<edge> & edges)
{
    set<edge>::iterator r = edges.find(ed);
    if(r != edges.end())
        return r->num;
    ed.num = edges.size();
    edges.insert(ed);
    return ed.num;
}

void VFEM::input_phys(const string & phys_filename)
{
    cout << "Reading physical parameters ..." << endl;

    // Чтение параметров физических областей
    ifstream phys_param;
    phys_param.open(phys_filename.c_str(), ios::in);
    if(!phys_param.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << phys_filename << endl;
        throw IO_FILE_ERROR;
    }
    double omega_global;
    phys_param >> omega_global;
    omega_global *= 2.0 * M_PI;

    size_t phys_num;
    phys_param >> phys_num;
    for(size_t i = 0; i < phys_num; i++)
    {
        phys_id id;
        phys_param >> id.gmsh_num >> id.type_of_element;
        phys_area * ph = &(phys[id]);

        ph->gmsh_num = id.gmsh_num;
        ph->type_of_elem = id.type_of_element;
        if(id.type_of_element == 4)
        {
            phys_param >> ph->mu >> ph->epsilon >> ph->sigma;
            ph->omega = omega_global;
            ph->type_of_bounds = 0;
            ph->J0 = 0.0;
        }
        else if(id.type_of_element == 2)
        {
            size_t parent_phys;
            phys_param >> parent_phys;
            phys_param >> ph->type_of_bounds;
            map<phys_id, phys_area>::const_iterator parent = phys.find(phys_id(4, parent_phys));
            if(parent != phys.end())
            {
                ph->omega = parent->second.omega;
                ph->mu = parent->second.mu;
                ph->epsilon = parent->second.epsilon;
                ph->sigma = parent->second.sigma;
                ph->J0 = parent->second.J0;
            }
            else
            {
                ph->omega = omega_global;
                ph->mu = 0.0;
                ph->epsilon = 0.0;
                ph->sigma = 0.0;
                ph->J0 = 0.0;
                cerr << "Warning: unaccounted parent \"" << parent_phys << "\" of phys area \""
                     << ph->gmsh_num << "\" (2), skipping..." << endl;
            }
        }
        else if(id.type_of_element == 1)
        {
            size_t parent_phys;
            phys_param >> parent_phys;
            phys_param >> ph->J0;
            map<phys_id, phys_area>::const_iterator parent = phys.find(phys_id(4, parent_phys));
            if(parent != phys.end())
            {
                ph->omega = parent->second.omega;
                ph->mu = parent->second.mu;
                ph->epsilon = parent->second.epsilon;
                ph->sigma = parent->second.sigma;
                ph->type_of_bounds = parent->second.type_of_bounds;
            }
            else
            {
                ph->omega = omega_global;
                ph->mu = 0.0;
                ph->epsilon = 0.0;
                ph->sigma = 0.0;
                ph->type_of_bounds = 0;
                cerr << "Warning: unaccounted parent \"" << parent_phys << "\" of phys area \""
                     << ph->gmsh_num << "\" (1), skipping..." << endl;
            }
        }
        else
        {
            cerr << "Warning: unaccounted type of element \"" << id.type_of_element
                 << "\" of phys area \"" << id.gmsh_num << "\", skipping..." << endl;
            string line;
            getline(phys_param, line);
        }
    }

    if(!phys_param.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << phys_filename << endl;
        throw IO_FILE_ERROR;
    }
    pss_num = 0;
    phys_param >> pss_num;
    if(pss_num)
        pss = new pair<point, cvector3> [pss_num];
    for(size_t i = 0; i < pss_num; i++)
    {
        for(size_t j = 0; j < 3; j++)
            phys_param >> pss[i].first[j];
        double d;
        for(size_t j = 0; j < 3; j++)
        {
            phys_param >> d;
            pss[i].second[j] = d;
        }
    }

    phys_param.close();
}

void VFEM::input_mesh(const string & gmsh_filename)
{
    cout << "Reading mesh ..." << endl;
    string line;

    // Чтение сетки
    ifstream gmsh_file;
    gmsh_file.open(gmsh_filename.c_str(), ios::in);

    cout << " > Reading nodes ..." << endl;

    do
        getline(gmsh_file, line);
    while(line.find("$Nodes") == string::npos && gmsh_file.good());

    if(!gmsh_file.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << gmsh_filename << endl;
        throw IO_FILE_ERROR;
    }

    // Чтение узлов
    gmsh_file >> nodes_num;
    nodes = new point [nodes_num];
    double max_coord[3] = {DBL_MIN, DBL_MIN, DBL_MIN};
    double min_coord[3] = {DBL_MAX, DBL_MAX, DBL_MAX};

    size_t fake_number;
    for(size_t i = 0; i < nodes_num; i++)
    {
        show_progress("", i, nodes_num);

        gmsh_file >> fake_number;
        for(size_t j = 0; j < 3; j++)
        {
            gmsh_file >> nodes[i][j];
            if(nodes[i][j] > max_coord[j])
                max_coord[j] = nodes[i][j];
            if(nodes[i][j] < min_coord[j])
                min_coord[j] = nodes[i][j];
        }
        nodes[i].num = i;
    }

    cout << " > Reading elements ..." << endl;

    do
        getline(gmsh_file, line);
    while(line.find("$Elements") == string::npos && gmsh_file.good());

    if(!gmsh_file.good())
    {
        cerr << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << gmsh_filename << endl;
        throw IO_FILE_ERROR;
    }

    vector<finite_element> fes_temp;
    vector<triangle> trs_temp;
    set<edge> edges_temp2;
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    set<edge> edges_surf_temp2;
#endif
    vector<edge_src> edges_src_temp;

    // Чтение конечных элементов
    gmsh_file >> fes_num;
    size_t type_of_elem;
    finite_element fake_element;
    triangle fake_triangle;
    edge_src fake_edge_src;
    vector<size_t> local_nodes_tet;
    local_nodes_tet.resize(4);
    vector<size_t> local_nodes_tr;
    local_nodes_tr.resize(3);

    for(size_t i = 0; i < fes_num; i++)
    {
        show_progress("", i, fes_num);

        gmsh_file >> fake_number >> type_of_elem;
        size_t phys_num, tags_num;
        gmsh_file >> tags_num;
        if(tags_num == 2)
            gmsh_file >> phys_num >> fake_number;
        else
        {
            cerr << "Warning: unaccounted number of tags, skipping..." << endl;
            for(size_t j = 0; j < tags_num; j++)
                gmsh_file >> fake_number;
        }

        if(type_of_elem == 4)
        {
            map<phys_id, phys_area>::iterator ph = phys.find(phys_id(4, phys_num));
            if(ph != phys.end())
                fake_element.phys = &(ph->second);
            else
            {
                cerr << "Error: can`t detect physical id " << phys_num << " in " << gmsh_filename << endl;
                throw IO_FILE_ERROR;
            }

            for(size_t j = 0; j < 4; j++)
                gmsh_file >> local_nodes_tet[j];
            sort(local_nodes_tet.begin(), local_nodes_tet.end());
            for(size_t j = 0; j < 4; j++)
                fake_element.nodes[j] = & nodes[local_nodes_tet[j] - 1];

            fake_element.edges[0] = (edge *) add_edge(edge(fake_element.nodes[0], fake_element.nodes[1]), edges_temp2);
            fake_element.edges[1] = (edge *) add_edge(edge(fake_element.nodes[0], fake_element.nodes[2]), edges_temp2);
            fake_element.edges[2] = (edge *) add_edge(edge(fake_element.nodes[0], fake_element.nodes[3]), edges_temp2);
            fake_element.edges[3] = (edge *) add_edge(edge(fake_element.nodes[1], fake_element.nodes[2]), edges_temp2);
            fake_element.edges[4] = (edge *) add_edge(edge(fake_element.nodes[1], fake_element.nodes[3]), edges_temp2);
            fake_element.edges[5] = (edge *) add_edge(edge(fake_element.nodes[2], fake_element.nodes[3]), edges_temp2);

            fes_temp.push_back(fake_element);
        }
        else if(type_of_elem == 2)
        {
            map<phys_id, phys_area>::iterator ph = phys.find(phys_id(2, phys_num));
            if(ph != phys.end())
                fake_triangle.phys = &(ph->second);
            else
            {
                cerr << "Error: can`t detect physical id " << phys_num << " in " << gmsh_filename << endl;
                throw IO_FILE_ERROR;
            }

            size_t bound_type = fake_triangle.get_phys_area().type_of_bounds;
            if(bound_type == 1)
                bound1_num++;
            else if(bound_type == 2)
                bound2_num++;
            else
            {
                cerr << "Error: unaccounted bound, breaking..." << endl;
                throw IO_FILE_ERROR;
            }

            for(size_t j = 0; j < 3; j++)
                gmsh_file >> local_nodes_tr[j];
            sort(local_nodes_tr.begin(), local_nodes_tr.end());
            for(size_t j = 0; j < 3; j++)
                fake_triangle.nodes[j] = & nodes[local_nodes_tr[j] - 1];

            fake_triangle.edges[0] = (edge *) add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[1]), edges_temp2);
            fake_triangle.edges[1] = (edge *) add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[2]), edges_temp2);
            fake_triangle.edges[2] = (edge *) add_edge(edge(fake_triangle.nodes[1], fake_triangle.nodes[2]), edges_temp2);
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
            if(bound_type == 1)
            {
                size_t e_num;
                e_num = add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[1]), edges_surf_temp2);
                fake_triangle.edges_surf[0] = e_num;
                e_num = add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[2]), edges_surf_temp2);
                fake_triangle.edges_surf[1] = e_num;
                e_num = add_edge(edge(fake_triangle.nodes[1], fake_triangle.nodes[2]), edges_surf_temp2);
                fake_triangle.edges_surf[2] = e_num;
            }
#endif
            trs_temp.push_back(fake_triangle);
        }
        else if(type_of_elem == 1)
        {
            map<phys_id, phys_area>::iterator ph = phys.find(phys_id(1, phys_num));
            if(ph != phys.end())
                fake_edge_src.phys = &(ph->second);
            else
            {
                cerr << "Error: can`t detect physical id " << phys_num << " in " << gmsh_filename << endl;
                throw IO_FILE_ERROR;
            }

            size_t nd[2];
            for(size_t j = 0; j < 2; j++)
                gmsh_file >> nd[j];

            if(nd[0] < nd[1])
            {
                fake_edge_src.direction = 1.0;
            }
            else
            {
                fake_edge_src.direction = -1.0;
                swap(nd[0], nd[1]);
            }

            for(size_t j = 0; j < 2; j++)
                fake_edge_src.nodes[j] = & nodes[nd[j] - 1];

            fake_edge_src.num = add_edge(edge(fake_edge_src.nodes[0], fake_edge_src.nodes[1]), edges_temp2);
            fake_edge_src.edge_main = NULL;
            edges_src_temp.push_back(fake_edge_src);
        }
        else
        {
            getline(gmsh_file, line);
        }
    }

    gmsh_file.close();

    cout << " > Converting data ..." << endl;

    edges_num = edges_temp2.size();
    edges = new edge [edges_num];
    size_t ied = 0;
    for(set<edge>::iterator i = edges_temp2.begin(); i != edges_temp2.end(); i++)
    {
        edges[i->num] = *i;
        show_progress("edges", ied++, edges_num);
    }

#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    if(bound1_num > 0)
    {
        edges_surf_num = edges_surf_temp2.size();
        for(set<edge>::iterator i = edges_surf_temp2.begin(); i != edges_surf_temp2.end(); i++)
        {
            set<edge>::iterator j = edges_temp2.find(*i);
            global_to_local[j->num] = i->num;
            global_to_local[j->num + edges_num] = i->num + edges_surf_num;
        }
        edges_surf_temp2.clear();
    }
#endif
    edges_temp2.clear();

    fes_num = fes_temp.size();
    if(fes_num > 0)
    {
        fes = new finite_element [fes_num];
        for(size_t i = 0; i < fes_num; i++)
        {
            show_progress("tetrahedrons", i, fes_num);
            fes[i] = fes_temp[i];
            for(size_t j = 0; j < 6; j++)
                fes[i].edges[j] = & edges[(size_t)fes[i].edges[j]];
            fes[i].init();
        }
        fes_temp.clear();
    }
    else
    {
        cerr << "Error: 0 tetrahedrons detected, breaking..." << endl;
        throw IO_FILE_ERROR;
    }

    edges_src_num = edges_src_temp.size();
    edges_src = new edge_src [edges_src_num];
    for(size_t i = 0; i < edges_src_num; i++)
    {
        show_progress("edges with source", i, edges_src_num);
        edges_src[i] = edges_src_temp[i];
        edges_src[i].edge_main = & edges[edges_src[i].num];
    }
    edges_src_temp.clear();

    trs_num = trs_temp.size();
    if(trs_num > 0)
    {
        trs = new triangle [trs_num];
        for(size_t i = 0; i < trs_num; i++)
        {
            show_progress("triangles", i, trs_num);
            trs[i] = trs_temp[i];
            for(size_t j = 0; j < 3; j++)
                trs[i].edges[j] = & edges[(size_t)trs[i].edges[j]];
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
            trs[i].init();
#else
            if(trs[i].get_phys_area().type_of_bounds == 1)
                for(size_t j = 0; j < 3; j++)
                {
                    edges_first.insert(trs[i].get_edge(j).num);
                    edges_first.insert(trs[i].get_edge(j).num + edges_num);
                }
#endif
        }
        trs_temp.clear();
    }
    else
    {
        cerr << "Error: 0 triangles detected, breaking..." << endl;
        throw IO_FILE_ERROR;
    }

    cout << " > Building tree ..." << endl;
    double diff_coord[3] =
    {
        max_coord[0] - min_coord[0],
        max_coord[1] - min_coord[1],
        max_coord[2] - min_coord[2]
    };
    for(size_t i = 0; i < 3; i++)
    {
        diff_coord[i] *= 1e-8;
        max_coord[i] += diff_coord[i];
        min_coord[i] -= diff_coord[i];
    }
    tree.make(min_coord[0], max_coord[0], min_coord[1], max_coord[1],
              min_coord[2], max_coord[2], fes, fes_num);

#if defined VFEM_USE_PML
    input_pml();
#endif
}

#if defined VFEM_USE_PML
void VFEM::input_pml()
{
    cout << " > Building PML-bound coordinates ..." << endl;
    // Границы не PML точек
    phys_pml.x0 = DBL_MAX;
    phys_pml.x1 = -DBL_MAX;
    phys_pml.y0 = DBL_MAX;
    phys_pml.y1 = -DBL_MAX;
    phys_pml.z0 = DBL_MAX;
    phys_pml.z1 = -DBL_MAX;
    map<point *, pair<cpoint *, finite_element *> > pml_nodes_tmp;

    for(size_t i = 0; i < fes_num; i++)
    {
        show_progress("scanned tetrahedrons", i, fes_num);
        if(is_pml(fes[i].barycenter, fes + i))
        {
            for(size_t j = 0; j < 4; j++)
                pml_nodes_tmp[fes[i].nodes[j]] = make_pair((cpoint *)NULL, fes + i);
        }
        else
        {
            for(size_t j = 0; j < 4; j++)
            {
                if(fes[i].nodes[j]->x > phys_pml.x1) phys_pml.x1 = fes[i].nodes[j]->x;
                if(fes[i].nodes[j]->y > phys_pml.y1) phys_pml.y1 = fes[i].nodes[j]->y;
                if(fes[i].nodes[j]->z > phys_pml.z1) phys_pml.z1 = fes[i].nodes[j]->z;
                if(fes[i].nodes[j]->x < phys_pml.x0) phys_pml.x0 = fes[i].nodes[j]->x;
                if(fes[i].nodes[j]->y < phys_pml.y0) phys_pml.y0 = fes[i].nodes[j]->y;
                if(fes[i].nodes[j]->z < phys_pml.z0) phys_pml.z0 = fes[i].nodes[j]->z;
            }
        }
    }
    cout << "  non-PML dimension: (" << phys_pml.x0 << "," << phys_pml.x1 << ")x(" << phys_pml.y0 << "," << phys_pml.y1 << ")x(" << phys_pml.z0 << "," << phys_pml.z1 << ")" << endl;

    // По умолчанию все координаты равны не-PML
    nodes_pml = new cpoint[nodes_num];
    for(size_t i = 0; i < nodes_num; i++)
        nodes_pml[i] = nodes[i];

    double pml_gauss_points_local[5] =
    {
        0.0,
        sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0,
        -sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0,
        sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0,
        -sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0,
    };
    double pml_gauss_weights[5] =
    {
        128.0 / 225.0,
        (322.0 + 13.0 * sqrt(70.0)) / 900.0,
        (322.0 + 13.0 * sqrt(70.0)) / 900.0,
        (322.0 - 13.0 * sqrt(70.0)) / 900.0,
        (322.0 - 13.0 * sqrt(70.0)) / 900.0
    };

    size_t i = 0;
    for(map<point *, pair<cpoint *, finite_element *> >::iterator it = pml_nodes_tmp.begin(); it != pml_nodes_tmp.end(); it++)
    {
        show_progress("replaced points", i++, pml_nodes_tmp.size());

        it->second.first = nodes_pml + (size_t)(it->first - nodes);
        double h[3] = {0, 0, 0}, beg[3] = {0, 0, 0};
        bool flag[3] = {false, false, false};
        point * p = it->first;
        finite_element * fefe = it->second.second;

        if(p->x > phys_pml.x1)
        {
            flag[0] = true;
            h[0] = p->x - phys_pml.x1;
            beg[0] = phys_pml.x1;
        }
        if(p->x < phys_pml.x0)
        {
            flag[0] = true;
            h[0] = p->x - phys_pml.x0;
            beg[0] = phys_pml.x0;
        }

        if(p->y > phys_pml.y1)
        {
            flag[1] = true;
            h[1] = p->y - phys_pml.y1;
            beg[1] = phys_pml.y1;
        }
        if(p->y < phys_pml.y0)
        {
            flag[1] = true;
            h[1] = p->y - phys_pml.y0;
            beg[1] = phys_pml.y0;
        }

        if(p->z > phys_pml.z1)
        {
            flag[2] = true;
            h[2] = p->z - phys_pml.z1;
            beg[2] = phys_pml.z1;
        }
        if(p->z < phys_pml.z0)
        {
            flag[2] = true;
            h[2] = p->z - phys_pml.z0;
            beg[2] = phys_pml.z0;
        }

        complex<double> new_point[3] = {0, 0, 0};

        size_t num_steps = 10;
        point h_small;
        for(size_t k = 0; k < 3; k++)
            h_small[k] = h[k] / (double)num_steps;

        for(size_t m = 0; m < num_steps; m++)
        {
            for(size_t k = 0; k < 5; k++)
            {
                point gauss_point_global((pml_gauss_points_local[k] + 1.0) / 2.0 * (h_small[0] * (double)m) + beg[0],
                                         (pml_gauss_points_local[k] + 1.0) / 2.0 * (h_small[1] * (double)m) + beg[1],
                                         (pml_gauss_points_local[k] + 1.0) / 2.0 * (h_small[2] * (double)m) + beg[2]);
                cvector3 s = get_s(gauss_point_global, fefe, & phys_pml);

                for(size_t j = 0; j < 3; j++)
                {
                    if(flag[j])
                        new_point[j] += s[j] * pml_gauss_weights[k];
                }
            }
        }

        for(size_t j = 0; j < 3; j++)
        {
            if(flag[j])
                (* (it->second.first))[j] = new_point[j] * h_small[j] / 2.0 + beg[j];
        }

        //cout << *(it->first) << " >> " << *(it->second.first) << endl;
    }

    for(size_t i = 0; i < fes_num; i++)
    {
        show_progress("re-init tetrahedrons", i, fes_num);
        for(size_t j = 0; j < 4; j++)
            fes[i].nodes_pml[j] = nodes_pml + (size_t)(fes[i].nodes[j] - nodes);
        fes[i].init_pml(get_s, & phys_pml);
    }
}
#endif
