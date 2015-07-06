#include "vfem.h"

size_t VFEM::add_edge(edge ed, set<edge> & edges_set)
{
    set<edge>::iterator r = edges_set.find(ed);
    if(r != edges_set.end())
        return r->num;
    ed.num = edges_set.size();
    edges_set.insert(ed);
    return ed.num;
}

#if BASIS_ORDER >= 2
size_t VFEM::add_face(face fc, set<face> & faces_set)
{
    set<face>::iterator r = faces_set.find(fc);
    if(r != faces_set.end())
        return r->num;
    fc.num = faces_set.size();
    faces_set.insert(fc);
    return fc.num;
}
#endif

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
            ph->mu *= consts::mu0;
            ph->epsilon *= consts::epsilon0;
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
                const phys_area * par = &(parent->second);
                ph->omega = par->omega;
                ph->mu = par->mu;
                ph->epsilon = par->epsilon;
                ph->sigma = par->sigma;
                ph->J0 = par->J0;
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
                const phys_area * par = &(parent->second);
                ph->omega = par->omega;
                ph->mu = par->mu;
                ph->epsilon = par->epsilon;
                ph->sigma = par->sigma;
                ph->type_of_bounds = par->type_of_bounds;
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
    size_t pss_num = 0;
    phys_param >> pss_num;
    if(pss_num)
        pss.resize(pss_num);
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
    size_t nodes_num;
    gmsh_file >> nodes_num;
    nodes.resize(nodes_num);
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

#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    set<edge> edges_surf_temp;
#if BASIS_ORDER >= 2
    set<face> faces_surf_temp;
#endif
#endif

    // Чтение конечных элементов
    size_t fes_num;
    gmsh_file >> fes_num;
    size_t type_of_elem = 0;
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
        size_t phys_num = 0, tags_num = 0;
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

            fake_element.edges[0] = (edge *) add_edge(edge(fake_element.nodes[0], fake_element.nodes[1]), edges);
            fake_element.edges[1] = (edge *) add_edge(edge(fake_element.nodes[0], fake_element.nodes[2]), edges);
            fake_element.edges[2] = (edge *) add_edge(edge(fake_element.nodes[0], fake_element.nodes[3]), edges);
            fake_element.edges[3] = (edge *) add_edge(edge(fake_element.nodes[1], fake_element.nodes[2]), edges);
            fake_element.edges[4] = (edge *) add_edge(edge(fake_element.nodes[1], fake_element.nodes[3]), edges);
            fake_element.edges[5] = (edge *) add_edge(edge(fake_element.nodes[2], fake_element.nodes[3]), edges);

#if BASIS_ORDER >= 2
            fake_element.faces[0] = (face *) add_face(face(fake_element.nodes[0], fake_element.nodes[1], fake_element.nodes[2]), faces);
            fake_element.faces[1] = (face *) add_face(face(fake_element.nodes[0], fake_element.nodes[1], fake_element.nodes[3]), faces);
            fake_element.faces[2] = (face *) add_face(face(fake_element.nodes[0], fake_element.nodes[2], fake_element.nodes[3]), faces);
            fake_element.faces[3] = (face *) add_face(face(fake_element.nodes[1], fake_element.nodes[2], fake_element.nodes[3]), faces);
#endif

            fes.push_back(fake_element);
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

            size_t bound_type = fake_triangle.phys->type_of_bounds;
            if(bound_type != 1 && bound_type != 2)
            {
                cerr << "Error: unaccounted bound, breaking..." << endl;
                throw IO_FILE_ERROR;
            }

            for(size_t j = 0; j < 3; j++)
                gmsh_file >> local_nodes_tr[j];
            sort(local_nodes_tr.begin(), local_nodes_tr.end());
            for(size_t j = 0; j < 3; j++)
                fake_triangle.nodes[j] = & nodes[local_nodes_tr[j] - 1];

            fake_triangle.edges[0] = (edge *) add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[1]), edges);
            fake_triangle.edges[1] = (edge *) add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[2]), edges);
            fake_triangle.edges[2] = (edge *) add_edge(edge(fake_triangle.nodes[1], fake_triangle.nodes[2]), edges);
#if BASIS_ORDER >= 2
            fake_triangle.faces = (face *) add_face(face(fake_triangle.nodes[0], fake_triangle.nodes[1], fake_triangle.nodes[2]), faces);
#endif
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
            if(bound_type == 1)
            {
                add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[1]), edges_surf_temp);
                add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[2]), edges_surf_temp);
                add_edge(edge(fake_triangle.nodes[1], fake_triangle.nodes[2]), edges_surf_temp);
#if BASIS_ORDER >= 2
                add_face(face(fake_triangle.nodes[0], fake_triangle.nodes[1], fake_triangle.nodes[2]), faces_surf_temp);
#endif
            }
#endif
            trs.push_back(fake_triangle);
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

            fake_edge_src.num = add_edge(edge(fake_edge_src.nodes[0], fake_edge_src.nodes[1]), edges);
            fake_edge_src.edge_main = NULL;
            edges_src.push_back(fake_edge_src);
        }
        else
        {
            getline(gmsh_file, line);
        }
    }

    gmsh_file.close();

#if  defined USE_CXX11
    /// WARNING: C++11
    fes.shrink_to_fit();
    trs.shrink_to_fit();
#endif

    cout << " > Converting data ..." << endl;

    // Индексируем ребра
    vector<edge *> edges_ind;
    edges_ind.resize(edges.size());
    size_t index_edge = 0;
    /// WARNING: const_cast для элемента set'а!
    for(set<edge>::iterator i = edges.begin(); i != edges.end(); ++i)
    {
        edge * ed = const_cast<edge *>(&(*i));
        edges_ind[ed->num] = ed;
        ed->num = index_edge++;
    }

#if BASIS_ORDER >= 2
    // Индексируем грани
    vector<face *> faces_ind;
    faces_ind.resize(faces.size());
    size_t index_face = 0;
    /// WARNING: const_cast для элемента set'а!
    for(set<face>::iterator i = faces.begin(); i != faces.end(); ++i)
    {
        face * fc = const_cast<face *>(&(*i));
        faces_ind[fc->num] = fc;
        fc->num = index_face++;
    }
#endif

    // Разбираемся с ребрами с источниками
    for(size_t i = 0; i < edges_src.size(); i++)
    {
        show_progress("edges with source", i, edges_src.size());
        edges_src[i].edge_main = edges_ind[edges_src[i].num];
        edges_src[i].num = edges_src[i].edge_main->num;
    }

    // Разбираемся с тетраэдрами
    if(fes.size() == 0)
    {
        cerr << "Error: 0 tetrahedrons detected, breaking..." << endl;
        throw IO_FILE_ERROR;
    }
    for(size_t i = 0; i < fes.size(); i++)
    {
        show_progress("tetrahedrons", i, fes.size());
        // Переносим ребра и грани
        for(size_t j = 0; j < 6; j++)
            fes[i].edges[j] = edges_ind[(size_t)fes[i].edges[j]];
#if BASIS_ORDER >= 2
        for(size_t j = 0; j < 4; j++)
            fes[i].faces[j] = faces_ind[(size_t)fes[i].faces[j]];
#endif
        // Заполняем степени свободы
        for(size_t j = 0; j < 6; j++)
            fes[i].dof[j] = fes[i].edges[j]->num;
        for(size_t j = 0; j < 4; j++)
            fes[i].ker_dof[j] = fes[i].nodes[j]->num;
#if BASIS_ORDER >= 2 || BASIS_TYPE == 2
        for(size_t j = 0; j < 6; j++)
            fes[i].dof[j + 6] = fes[i].edges[j]->num + edges.size();
        for(size_t j = 0; j < 6; j++)
            fes[i].ker_dof[j + 4] = fes[i].edges[j]->num + nodes.size();
#endif
#if BASIS_ORDER >= 2
        for(size_t j = 0; j < 4; j++)
            fes[i].dof[j + 12] = fes[i].faces[j]->num + 2 * edges.size();
        for(size_t j = 0; j < 4; j++)
            fes[i].dof[j + 16] = fes[i].faces[j]->num + 2 * edges.size() + faces.size();
#endif
#if BASIS_ORDER > 2 || (BASIS_TYPE == 2 && BASIS_ORDER == 2)
        for(size_t j = 0; j < 4; j++)
            fes[i].dof[j + 20] = fes[i].faces[j]->num + 2 * edges.size() + 2 * faces.size();
        for(size_t j = 0; j < 6; j++)
            fes[i].dof[j + 24] = fes[i].edges[j]->num + 2 * edges.size() + 3 * faces.size();
        for(size_t j = 0; j < 4; j++)
            fes[i].ker_dof[j + 10] = fes[i].faces[j]->num + edges.size() + nodes.size();
        for(size_t j = 0; j < 6; j++)
            fes[i].ker_dof[j + 14] = fes[i].edges[j]->num + faces.size() + edges.size() + nodes.size();;
#endif
        // Инициализируем
        fes[i].init();
    }

    // Разбираемся с треугольниками
    if(trs.size() == 0)
    {
        cerr << "Error: 0 triangles detected, breaking..." << endl;
        throw IO_FILE_ERROR;
    }
    for(size_t i = 0; i < trs.size(); i++)
    {
        show_progress("triangles", i, trs.size());
        for(size_t j = 0; j < 3; j++)
            trs[i].edges[j] = edges_ind[(size_t)trs[i].edges[j]];
#if BASIS_ORDER >= 2
        trs[i].faces = faces_ind[(size_t)trs[i].faces];
#endif

        // Заполняем степени свободы
        // Первый неполный
        for(size_t j = 0; j < 3; j++)
            trs[i].dof[j] = trs[i].edges[j]->num;
        for(size_t j = 0; j < 3; j++)
            trs[i].ker_dof[j] = trs[i].nodes[j]->num;
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
        if(trs[i].phys->type_of_bounds == 1)
            for(size_t j = 0; j < 3; j++)
                trs[i].dof_surf[j] = edges_surf_temp.find(* trs[i].edges[j])->num;
#endif
        // Первый полный
#if BASIS_ORDER >= 2 || BASIS_TYPE == 2
        for(size_t j = 0; j < 3; j++)
            trs[i].dof[j + 3] = trs[i].edges[j]->num + edges.size();
        for(size_t j = 0; j < 3; j++)
            trs[i].ker_dof[j + 3] = trs[i].edges[j]->num + nodes.size();
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
        if(trs[i].phys->type_of_bounds == 1)
            for(size_t j = 0; j < 3; j++)
                trs[i].dof_surf[j + 3] = edges_surf_temp.find(* trs[i].edges[j])->num + edges_surf_temp.size();
#endif
#endif
        // Второй неполный
#if BASIS_ORDER >= 2
        trs[i].dof[6] = trs[i].faces->num + 2 * edges.size();
        trs[i].dof[7] = trs[i].faces->num + 2 * edges.size() + faces.size();
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
        if(trs[i].phys->type_of_bounds == 1)
        {
            trs[i].dof_surf[6] = faces_surf_temp.find(* trs[i].faces)->num + 2 * edges_surf_temp.size();
            trs[i].dof_surf[7] = faces_surf_temp.find(* trs[i].faces)->num + 2 * edges_surf_temp.size() + faces_surf_temp.size();
        }
#endif
#endif
        // Второй полный
#if BASIS_ORDER > 2 || (BASIS_TYPE == 2 && BASIS_ORDER == 2)
        trs[i].dof[8] = trs[i].faces->num + 2 * edges.size() + 2 * faces.size();
        for(size_t j = 0; j < 3; j++)
            trs[i].dof[j + 9] = trs[i].edges[j]->num + 2 * edges.size() + 3 * faces.size();
        trs[i].ker_dof[6] = trs[i].faces->num + edges.size() + nodes.size();
        for(size_t j = 0; j < 3; j++)
            trs[i].ker_dof[j + 7] = trs[i].edges[j]->num + edges.size() + faces.size() + nodes.size();
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
        if(trs[i].phys->type_of_bounds == 1)
        {
            trs[i].dof_surf[8] = faces_surf_temp.find(* trs[i].faces)->num + 2 * edges_surf_temp.size() + 2 * faces_surf_temp.size();
            for(size_t j = 0; j < 3; j++)
                trs[i].dof_surf[j + 9] = edges_surf_temp.find(* trs[i].edges[j])->num + 2 * edges_surf_temp.size() + 3 * faces_surf_temp.size();
        }
#endif
#endif

        // Инициализируем
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
        if(trs[i].phys->type_of_bounds == 1)
            for(size_t j = 0; j < basis::tr_bf_num; j++)
                global_to_local[trs[i].dof[j]] = trs[i].dof_surf[j];
        trs[i].init();
#else
        if(trs[i].phys->type_of_bounds == 1)
            for(size_t j = 0; j < basis::tr_bf_num; j++)
                dof_first.insert(trs[i].dof[j]);
#endif
        if(trs[i].phys->type_of_bounds == 1)
            for(size_t j = 0; j < basis::tr_ker_bf_num; j++)
                ker_dof_first.insert(trs[i].ker_dof[j]);
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
#if defined USE_CXX11
    /// WARNING: C++11
    tree.make(min_coord[0], max_coord[0], min_coord[1], max_coord[1],
              min_coord[2], max_coord[2], fes.data(), fes.size());
#else
    /// WARNING: &(fes[0]) для вектора!
    tree.make(min_coord[0], max_coord[0], min_coord[1], max_coord[1],
              min_coord[2], max_coord[2], &(fes[0]), fes.size());
#endif

#if defined VFEM_USE_PML
    input_pml();
#endif

#if BASIS_ORDER == 1 && BASIS_TYPE == 1
    dof_num = edges.size();
    ker_dof_num = nodes.size();
#endif
#if BASIS_ORDER == 1 && BASIS_TYPE == 2
    dof_num = 2 * edges.size();
    ker_dof_num = nodes.size() + edges.size();
#endif
#if BASIS_ORDER == 2 && BASIS_TYPE == 1
    dof_num = 2 * edges.size() + 2 * faces.size();
    ker_dof_num = nodes.size() + edges.size();
#endif
#if BASIS_ORDER == 2 && BASIS_TYPE == 2
    dof_num = 3 * edges.size() + 3 * faces.size();
    ker_dof_num = nodes.size() + 2 * edges.size() + faces.size();
#endif

    cout << "Statistics:" << endl;
    cout << " # Tetrehedrons: " << fes.size() << endl;
    cout << " # Nodes:        " << nodes.size() << endl;
    cout << " # Edges:        " << edges.size() << endl;
#if BASIS_ORDER >= 2
    cout << " # Faces:        " << faces.size() << endl;
#endif
    cout << " # SLAE size:    " << dof_num << endl;
    cout << " # SLAE ker:     " << ker_dof_num << endl;
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    cout << " # SLAE surf:    " << global_to_local.size() << endl;
#endif
}

#if defined VFEM_USE_PML
cpoint VFEM::convert_point_to_pml(const point * p, const finite_element * fefe) const
{
    static const double pml_gauss_points_local[5] =
    {
        0.0,
        sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0,
        -sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0,
        sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0,
        -sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0,
    };
    static const double pml_gauss_weights[5] =
    {
        128.0 / 225.0,
        (322.0 + 13.0 * sqrt(70.0)) / 900.0,
        (322.0 + 13.0 * sqrt(70.0)) / 900.0,
        (322.0 - 13.0 * sqrt(70.0)) / 900.0,
        (322.0 - 13.0 * sqrt(70.0)) / 900.0
    };

    double h[3] = {0, 0, 0}, beg[3] = {0, 0, 0};
    bool flag[3] = {false, false, false};

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

    cpoint cp(* p);
    for(size_t j = 0; j < 3; j++)
    {
        if(flag[j])
            cp[j] = new_point[j] * h_small[j] / 2.0 + beg[j];
    }
    return cp;
}

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
    map<point *, pair<cpoint, finite_element *> > pml_nodes_cache;

    for(size_t i = 0; i < fes.size(); i++)
    {
        show_progress("scanned tetrahedrons", i, fes.size());
        finite_element * fes_i = &(fes[i]);
        if(is_pml(fes_i->barycenter, fes_i))
        {
            for(size_t j = 0; j < 4; j++)
                pml_nodes_cache[fes_i->nodes[j]] = make_pair(cpoint(), fes_i);
        }
        else
        {
            for(size_t j = 0; j < 4; j++)
            {
                point * nodes_j = fes_i->nodes[j];
                if(nodes_j->x > phys_pml.x1) phys_pml.x1 = nodes_j->x;
                if(nodes_j->y > phys_pml.y1) phys_pml.y1 = nodes_j->y;
                if(nodes_j->z > phys_pml.z1) phys_pml.z1 = nodes_j->z;
                if(nodes_j->x < phys_pml.x0) phys_pml.x0 = nodes_j->x;
                if(nodes_j->y < phys_pml.y0) phys_pml.y0 = nodes_j->y;
                if(nodes_j->z < phys_pml.z0) phys_pml.z0 = nodes_j->z;
            }
        }
    }
    cout << "  non-PML dimension: (" << phys_pml.x0 << "," << phys_pml.x1 << ")x(" << phys_pml.y0 << "," << phys_pml.y1 << ")x(" << phys_pml.z0 << "," << phys_pml.z1 << ")" << endl;

    size_t i = 0;
    for(map<point *, pair<cpoint, finite_element *> >::iterator it = pml_nodes_cache.begin(); it != pml_nodes_cache.end(); ++it)
    {
        show_progress("replaced points", i++, pml_nodes_cache.size());

        point * p = it->first;
        finite_element * fefe = it->second.second;
        it->second.first = convert_point_to_pml(p, fefe);
    }

    for(size_t i = 0; i < fes.size(); i++)
    {
        show_progress("re-init tetrahedrons", i, fes.size());
        cpoint cp[4];
        size_t ph_curr = fes[i].phys->gmsh_num;
        for(size_t j = 0; j < 4; j++)
        {
            if(is_pml(fes[i].barycenter, &(fes[i])))
            {
                map<point *, pair<cpoint, finite_element *> >::iterator it = pml_nodes_cache.find(fes[i].nodes[j]);
                // Если есть в кэше, то возьмем из него
                if(it->second.second->phys->gmsh_num == ph_curr)
                    cp[j] = it->second.first;
                // А иначе рассчитаем
                else
                    cp[j] = convert_point_to_pml(fes[i].nodes[j], &(fes[i]));
            }
            else
            {
                cp[j] = cpoint(* (fes[i].nodes[j]));
            }
        }
        fes[i].init_pml(get_s, & phys_pml, cp);
    }
}
#endif
