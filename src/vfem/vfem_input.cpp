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
            }
            else
            {
                ph->omega = omega_global;
                ph->mu = 0.0;
                ph->epsilon = 0.0;
                ph->sigma = 0.0;
                cerr << "Warning: unaccounted parent \"" << parent_phys << "\" of phys area \""
                     << ph->gmsh_num << "\" (2), skipping..." << endl;
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

    set<edge> edges_surf_temp;
#if BASIS_ORDER >= 2
    set<face> faces_surf_temp;
#endif

    // Чтение конечных элементов
    size_t fes_num;
    gmsh_file >> fes_num;
    size_t type_of_elem = 0;
    finite_element fake_element;
    triangle fake_triangle;
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
            if(bound_type == 1)
            {
                add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[1]), edges_surf_temp);
                add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[2]), edges_surf_temp);
                add_edge(edge(fake_triangle.nodes[1], fake_triangle.nodes[2]), edges_surf_temp);
#if BASIS_ORDER >= 2
                add_face(face(fake_triangle.nodes[0], fake_triangle.nodes[1], fake_triangle.nodes[2]), faces_surf_temp);
#endif
            }
            trs.push_back(fake_triangle);
        }
        else
        {
            getline(gmsh_file, line);
        }
    }
    gmsh_file.close();

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
#if BASIS_ORDER >= 2 || BASIS_TYPE == 2
        for(size_t j = 0; j < 6; j++)
            fes[i].dof[j + 6] = fes[i].edges[j]->num + edges.size();
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
        if(trs[i].phys->type_of_bounds == 1)
            for(size_t j = 0; j < 3; j++)
                trs[i].dof_surf[j] = edges_surf_temp.find(* trs[i].edges[j])->num;
        // Первый полный
#if BASIS_ORDER >= 2 || BASIS_TYPE == 2
        for(size_t j = 0; j < 3; j++)
            trs[i].dof[j + 3] = trs[i].edges[j]->num + edges.size();
        if(trs[i].phys->type_of_bounds == 1)
            for(size_t j = 0; j < 3; j++)
                trs[i].dof_surf[j + 3] = edges_surf_temp.find(* trs[i].edges[j])->num + edges_surf_temp.size();
#endif
        // Второй неполный
#if BASIS_ORDER >= 2
        trs[i].dof[6] = trs[i].faces->num + 2 * edges.size();
        trs[i].dof[7] = trs[i].faces->num + 2 * edges.size() + faces.size();
        if(trs[i].phys->type_of_bounds == 1)
        {
            trs[i].dof_surf[6] = faces_surf_temp.find(* trs[i].faces)->num + 2 * edges_surf_temp.size();
            trs[i].dof_surf[7] = faces_surf_temp.find(* trs[i].faces)->num + 2 * edges_surf_temp.size() + faces_surf_temp.size();
        }
#endif
        // Второй полный
#if BASIS_ORDER > 2 || (BASIS_TYPE == 2 && BASIS_ORDER == 2)
        trs[i].dof[8] = trs[i].faces->num + 2 * edges.size() + 2 * faces.size();
        for(size_t j = 0; j < 3; j++)
            trs[i].dof[j + 9] = trs[i].edges[j]->num + 2 * edges.size() + 3 * faces.size();
        if(trs[i].phys->type_of_bounds == 1)
        {
            trs[i].dof_surf[8] = faces_surf_temp.find(* trs[i].faces)->num + 2 * edges_surf_temp.size() + 2 * faces_surf_temp.size();
            for(size_t j = 0; j < 3; j++)
                trs[i].dof_surf[j + 9] = edges_surf_temp.find(* trs[i].edges[j])->num + 2 * edges_surf_temp.size() + 3 * faces_surf_temp.size();
        }
#endif

        // Инициализируем
        if(trs[i].phys->type_of_bounds == 1)
            for(size_t j = 0; j < basis::tr_bf_num; j++)
                global_to_local[trs[i].dof[j]] = trs[i].dof_surf[j];
        trs[i].init();
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
    /// WARNING: &(fes[0]) для вектора!
    tree.make(min_coord[0], max_coord[0], min_coord[1], max_coord[1],
              min_coord[2], max_coord[2], &(fes[0]), fes.size());

#if BASIS_ORDER == 1 && BASIS_TYPE == 1
    dof_num = edges.size();
#endif
#if BASIS_ORDER == 1 && BASIS_TYPE == 2
    dof_num = 2 * edges.size();
#endif
#if BASIS_ORDER == 2 && BASIS_TYPE == 1
    dof_num = 2 * edges.size() + 2 * faces.size();
#endif
#if BASIS_ORDER == 2 && BASIS_TYPE == 2
    dof_num = 3 * edges.size() + 3 * faces.size();
#endif

    cout << "Statistics:" << endl;
    cout << " # Tetrehedrons: " << fes.size() << endl;
    cout << " # Nodes:        " << nodes.size() << endl;
    cout << " # Edges:        " << edges.size() << endl;
#if BASIS_ORDER >= 2
    cout << " # Faces:        " << faces.size() << endl;
#endif
    cout << " # SLAE size:    " << dof_num << endl;
    cout << " # SLAE surf:    " << global_to_local.size() << endl;
}
