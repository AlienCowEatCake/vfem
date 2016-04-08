#include "vfem.h"

// gmsh-2.11.0-source/Common/GmshDefines.h
#define MSH_LIN_2    1  // 2-node line.
#define MSH_TRI_3    2  // 3-node triangle.
#define MSH_TET_4    4  // 4-node tetrahedron.

template<typename U, typename V>
void push_back_wrapper(vector<U> & v, const V & e)
{
    size_t size = v.size();
    size_t capacity = v.capacity();
    size_t danger = 268435456 / sizeof(U); // 256 MiB
    // Грязный хак
    if(capacity == size)
    {
        if(size > danger)
            v.reserve((size_t)((double)capacity * 1.5));
        else
            v.reserve(capacity * 2);
    }
    // А теперь уже и добавим
    v.push_back(e);
}

size_t VFEM::add_edge(edge ed, set<edge> & edges_set)
{
    set<edge>::iterator r = edges_set.find(ed);
    if(r != edges_set.end())
        return r->num;
    ed.num = edges_set.size();
    edges_set.insert(ed);
    return ed.num;
}

size_t VFEM::add_face(face fc, set<face> & faces_set)
{
    set<face>::iterator r = faces_set.find(fc);
    if(r != faces_set.end())
        return r->num;
    fc.num = faces_set.size();
    faces_set.insert(fc);
    return fc.num;
}

bool VFEM::input_phys(const string & phys_filename)
{
    cout << "Reading physical parameters ..." << endl;

    // Чтение параметров физических областей
    ifstream phys_param;
    phys_param.open(phys_filename.c_str(), ios::in);
    if(!phys_param.good())
    {
        cout << "[Phys Config] Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << phys_filename << endl;
        return false;
    }

    double mu_default    = consts::mu0;
    double eps_default   = consts::epsilon0;
    double sigma_default = 0.0;
    double omega_global  = -1.0 * 2.0 * M_PI;

    string to_lowercase(const string & str);
    string trim(const string & str);

    string line;
    getline(phys_param, line);
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

            if(section == "global")
            {
                do
                {
                    getline(phys_param, line);
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
                            if(param == "frequency")  omega_global = tmp * 2.0 * M_PI;
                            else if(param == "mu")    mu_default = tmp * consts::mu0;
                            else if(param == "eps")   eps_default = tmp * consts::epsilon0;
                            else if(param == "sigma") sigma_default = tmp;
                            else cout << "[Phys Config] Unsupported param \"" << param << "\" in section \""
                                      << section << (subsection.empty() ? string("") : (string(".") + subsection))
                                      << "\"" << endl;
                            //cout << "  param = " << param << endl;
                            //cout << "  value = " << value << endl;
                        }
                    }
                }
                while(phys_param.good() && !(line.length() > 1 && line[0] == '['));
            }
            else if(section == "phys3d")
            {
                if(!subsection.empty())
                {
                    if(omega_global < 0)
                    {
                        cout << "[Phys Config] Unknown frequency in section \"" << section
                             << "." << subsection << "\"" << endl;
                        return false;
                    }

                    stringstream sst(subsection);
                    size_t gmsh_num;
                    sst >> gmsh_num;
                    phys_area * ph = &(phys[phys_id(MSH_TET_4, gmsh_num)]);

                    ph->gmsh_num = gmsh_num;
                    ph->type_of_elem = MSH_TET_4;
                    ph->mu = mu_default;
                    ph->epsilon = eps_default;
                    ph->sigma = sigma_default;
                    ph->omega = omega_global;
                    ph->type_of_bounds = 0;
                    ph->J0 = 0.0;
                    ph->E0 = 0.0;

                    do
                    {
                        getline(phys_param, line);
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
                                if(param == "mu")         ph->mu = tmp * consts::mu0;
                                else if(param == "eps")   ph->epsilon = tmp * consts::epsilon0;
                                else if(param == "sigma") ph->sigma = tmp;
                                else cout << "[Phys Config] Unsupported param \"" << param << "\" in section \""
                                          << section << (subsection.empty() ? string("") : (string(".") + subsection))
                                          << "\"" << endl;
                                //cout << "  param = " << param << endl;
                                //cout << "  value = " << value << endl;
                            }
                        }
                    }
                    while(phys_param.good() && !(line.length() > 1 && line[0] == '['));
                }
                else
                {
                    cout << "[Phys Config] No subsection in section \"" << section << "\"" << endl;
                    return false;
                }
            }
            else if(section == "phys2d")
            {
                if(!subsection.empty())
                {
                    if(omega_global < 0)
                    {
                        cout << "[Phys Config] Unknown frequency in section \"" << section
                             << "." << subsection << "\"" << endl;
                        return false;
                    }

                    stringstream sst(subsection);
                    size_t gmsh_num;
                    sst >> gmsh_num;
                    phys_area * ph = &(phys[phys_id(MSH_TRI_3, gmsh_num)]);

                    ph->gmsh_num = gmsh_num;
                    ph->type_of_elem = MSH_TRI_3;
                    ph->mu = mu_default;
                    ph->epsilon = eps_default;
                    ph->sigma = sigma_default;
                    ph->omega = omega_global;
                    ph->type_of_bounds = 2;
                    ph->J0 = 0.0;
                    ph->E0 = 0.0;

                    do
                    {
                        getline(phys_param, line);
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
                                size_t tmp;
                                sst >> tmp;
                                if(param == "boundary") ph->type_of_bounds = tmp;
                                else if(param == "parent")
                                {
                                    map<phys_id, phys_area>::const_iterator parent =
                                            phys.find(phys_id(MSH_TET_4, tmp));
                                    if(parent != phys.end())
                                    {
                                        const phys_area * par = &(parent->second);
                                        ph->mu = par->mu;
                                        ph->epsilon = par->epsilon;
                                        ph->sigma = par->sigma;
                                    }
                                    else
                                    {
                                        cout << "[Phys Config] Warning: unaccounted parent \""
                                             << tmp << "\" of phys area \"" << ph->gmsh_num
                                             << "\" (2 - MSH_TRI_3), skipping..." << endl;
                                    }
                                }
                                else cout << "[Phys Config] Unsupported param \"" << param << "\" in section \""
                                          << section << (subsection.empty() ? string("") : (string(".") + subsection))
                                          << "\"" << endl;
                                //cout << "  param = " << param << endl;
                                //cout << "  value = " << value << endl;
                            }
                        }
                    }
                    while(phys_param.good() && !(line.length() > 1 && line[0] == '['));
                }
                else
                {
                    cout << "[Phys Config] No subsection in section \"" << section << "\"" << endl;
                    return false;
                }
            }
            else if(section == "phys1d")
            {
                if(!subsection.empty())
                {
                    if(omega_global < 0)
                    {
                        cout << "[Phys Config] Unknown frequency in section \"" << section
                             << "." << subsection << "\"" << endl;
                        return false;
                    }

                    stringstream sst(subsection);
                    size_t gmsh_num;
                    sst >> gmsh_num;
                    phys_area * ph = &(phys[phys_id(MSH_LIN_2, gmsh_num)]);

                    ph->gmsh_num = gmsh_num;
                    ph->type_of_elem = MSH_LIN_2;
                    ph->mu = mu_default;
                    ph->epsilon = eps_default;
                    ph->sigma = sigma_default;
                    ph->omega = omega_global;
                    ph->type_of_bounds = 0;
                    ph->J0 = 0.0;
                    ph->E0 = 0.0;

                    do
                    {
                        getline(phys_param, line);
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
                                if(param == "current") sst >> ph->J0;
                                else if(param == "parent")
                                {
                                    size_t parent_phys;
                                    sst >> parent_phys;
                                    map<phys_id, phys_area>::const_iterator parent =
                                            phys.find(phys_id(MSH_TET_4, parent_phys));
                                    if(parent != phys.end())
                                    {
                                        const phys_area * par = &(parent->second);
                                        ph->mu = par->mu;
                                        ph->epsilon = par->epsilon;
                                        ph->sigma = par->sigma;
                                    }
                                    else
                                    {
                                        cout << "[Phys Config] Warning: unaccounted parent \""
                                             << parent_phys << "\" of phys area \"" << ph->gmsh_num
                                             << "\" (1 - MSH_LIN_2), skipping..." << endl;
                                    }
                                }
                                else cout << "[Phys Config] Unsupported param \"" << param << "\" in section \""
                                          << section << (subsection.empty() ? string("") : (string(".") + subsection))
                                          << "\"" << endl;
                                //cout << "  param = " << param << endl;
                                //cout << "  value = " << value << endl;
                            }
                        }
                    }
                    while(phys_param.good() && !(line.length() > 1 && line[0] == '['));
                }
                else
                {
                    cout << "[Phys Config] No subsection in section \"" << section << "\"" << endl;
                    return false;
                }
            }
            else if(section == "phys0d")
            {
                point pss_point;
                cvector3 pss_value;
                do
                {
                    getline(phys_param, line);
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
                            if(param == "point_x")      pss_point.x = tmp;
                            else if(param == "point_y") pss_point.y = tmp;
                            else if(param == "point_z") pss_point.z = tmp;
                            else if(param == "value_x") pss_value.x = tmp;
                            else if(param == "value_y") pss_value.y = tmp;
                            else if(param == "value_z") pss_value.z = tmp;
                            else cout << "[Phys Config] Unsupported param \"" << param << "\" in section \""
                                      << section << (subsection.empty() ? string("") : (string(".") + subsection))
                                      << "\"" << endl;
                            //cout << "  param = " << param << endl;
                            //cout << "  value = " << value << endl;
                        }
                    }
                }
                while(phys_param.good() && !(line.length() > 1 && line[0] == '['));
                push_back_wrapper(pss, make_pair(pss_point, pss_value));
            }
            else if(section == "electrode")
            {
                if(!subsection.empty())
                {
                    if(omega_global < 0)
                    {
                        cout << "[Phys Config] Unknown frequency in section \"" << section
                             << "." << subsection << "\"" << endl;
                        return false;
                    }

                    stringstream sst(subsection);
                    size_t gmsh_num;
                    sst >> gmsh_num;
                    phys_area * ph = &(phys[phys_id(MSH_LIN_2, gmsh_num)]);

                    ph->gmsh_num = gmsh_num;
                    ph->type_of_elem = MSH_LIN_2;
                    ph->mu = mu_default;
                    ph->epsilon = eps_default;
                    ph->sigma = sigma_default;
                    ph->omega = omega_global;
                    ph->type_of_bounds = 1;
                    ph->J0 = 0.0;
                    ph->E0 = 0.0;

                    double I = 1.0;
                    double l = 1.0;
                    double rho = 1.0;
                    double S = 1.0;

                    do
                    {
                        getline(phys_param, line);
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
                                if     (param == "i")   sst >> I;
                                else if(param == "l")   sst >> l;
                                else if(param == "rho") sst >> rho;
                                else if(param == "s")   sst >> S;
                                else if(param == "parent")
                                {
                                    size_t parent_phys;
                                    sst >> parent_phys;
                                    map<phys_id, phys_area>::const_iterator parent =
                                            phys.find(phys_id(MSH_TET_4, parent_phys));
                                    if(parent != phys.end())
                                    {
                                        const phys_area * par = &(parent->second);
                                        ph->mu = par->mu;
                                        ph->epsilon = par->epsilon;
                                        ph->sigma = par->sigma;
                                    }
                                    else
                                    {
                                        cout << "[Phys Config] Warning: unaccounted parent \""
                                             << parent_phys << "\" of phys area \"" << ph->gmsh_num
                                             << "\" (1 - MSH_LIN_2), skipping..." << endl;
                                    }
                                }
                                else cout << "[Phys Config] Unsupported param \"" << param << "\" in section \""
                                          << section << (subsection.empty() ? string("") : (string(".") + subsection))
                                          << "\"" << endl;
                                //cout << "  param = " << param << endl;
                                //cout << "  value = " << value << endl;
                            }
                        }
                    }
                    while(phys_param.good() && !(line.length() > 1 && line[0] == '['));

                    double R = rho * l / S;
                    double U = I * R;
                    ph->E0 = U / l;
                }
                else
                {
                    cout << "[Phys Config] No subsection in section \"" << section << "\"" << endl;
                    return false;
                }
            }
        }
        else
        {
            getline(phys_param, line);
            line = trim(line);
        }
    }
    while(phys_param.good());

    phys_param.close();


    // Заполним данные и в вычисляемых функциях
    for(map<size_t, array_t<evaluator_helmholtz, 3> >::iterator
        it = config.analytical.values.begin(); it != config.analytical.values.end(); ++it)
    {
        phys_area * ph = &(phys.find(phys_id(MSH_TET_4, it->first))->second);
        complex<double> k2(- ph->omega * ph->omega * ph->epsilon,
                           ph->omega * ph->sigma);
        for(size_t i = 0; i < 3; i++)
        {
            evaluator_helmholtz * ev_curr = &(it->second[i]);
            ev_curr->set_J0(ph->J0);
            ev_curr->set_omega(omega_global);
            ev_curr->set_epsilon(ph->epsilon);
            ev_curr->set_sigma(ph->sigma);
            ev_curr->set_mu(ph->mu);
            ev_curr->set_k2(k2);
            ev_curr->simplify();
            switch(config.jit_type)
            {
            case evaluator3::JIT_EXTCALL:
                ev_curr->compile_extcall();
                break;
            case evaluator3::JIT_INLINE:
                ev_curr->compile_inline();
                break;
            default:
                break;
            }
        }
    }
    for(map<size_t, array_t<evaluator_helmholtz, 3> >::iterator
        it = config.boundary.values.begin(); it != config.boundary.values.end(); ++it)
    {
        phys_area * ph = &(phys.find(phys_id(MSH_TRI_3, it->first))->second);
        complex<double> k2(- ph->omega * ph->omega * ph->epsilon,
                           ph->omega * ph->sigma);
        for(size_t i = 0; i < 3; i++)
        {
            evaluator_helmholtz * ev_curr = &(it->second[i]);
            ev_curr->set_J0(ph->J0);
            ev_curr->set_omega(omega_global);
            ev_curr->set_epsilon(ph->epsilon);
            ev_curr->set_sigma(ph->sigma);
            ev_curr->set_mu(ph->mu);
            ev_curr->set_k2(k2);
            ev_curr->simplify();
            switch(config.jit_type)
            {
            case evaluator3::JIT_EXTCALL:
                ev_curr->compile_extcall();
                break;
            case evaluator3::JIT_INLINE:
                ev_curr->compile_inline();
                break;
            default:
                break;
            }
        }
    }
    for(map<size_t, array_t<evaluator_helmholtz, 3> >::iterator
        it = config.right.values.begin(); it != config.right.values.end(); ++it)
    {
        phys_area * ph = &(phys.find(phys_id(MSH_TET_4, it->first))->second);
        complex<double> k2(- ph->omega * ph->omega * ph->epsilon,
                           ph->omega * ph->sigma);
        for(size_t i = 0; i < 3; i++)
        {
            evaluator_helmholtz * ev_curr = &(it->second[i]);
            ev_curr->set_J0(ph->J0);
            ev_curr->set_omega(omega_global);
            ev_curr->set_epsilon(ph->epsilon);
            ev_curr->set_sigma(ph->sigma);
            ev_curr->set_mu(ph->mu);
            ev_curr->set_k2(k2);
            ev_curr->simplify();
            switch(config.jit_type)
            {
            case evaluator3::JIT_EXTCALL:
                ev_curr->compile_extcall();
                break;
            case evaluator3::JIT_INLINE:
                ev_curr->compile_inline();
                break;
            default:
                break;
            }
        }
    }
    for(size_t i = 0; i < 3; i++)
    {
        config.analytical.default_value[i].set_omega(omega_global);
        config.boundary.default_value[i].set_omega(omega_global);
        config.right.default_value[i].set_omega(omega_global);
    }

    return true;
}

bool VFEM::input_mesh(const string & gmsh_filename)
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
        cout << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << gmsh_filename << endl;
        return false;
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
        cout << "Error in " << __FILE__ << ":" << __LINE__
             << " while reading file " << gmsh_filename << endl;
        return false;
    }

    set<edge> edges_surf_temp;
    set<face> faces_surf_temp;

    // Чтение конечных элементов
    size_t fes_num;
    gmsh_file >> fes_num;
    size_t type_of_elem = 0;
    finite_element fake_element;
    triangle_base fake_triangle;
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
            cout << "Warning: unaccounted number of tags, skipping..." << endl;
            for(size_t j = 0; j < tags_num; j++)
                gmsh_file >> fake_number;
        }

        if(type_of_elem == MSH_TET_4)
        {
            map<phys_id, phys_area>::iterator ph = phys.find(phys_id(MSH_TET_4, phys_num));
            if(ph != phys.end())
                fake_element.phys = &(ph->second);
            else
            {
                cout << "Error in " << __FILE__ << ":" << __LINE__
                     << " while reading file " << gmsh_filename << endl;
                cout << "Can`t detect physical id " << phys_num << " in 4 - MSH_TET_4 (tetrahedron)" << endl;
                return false;
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

            if(config.basis.order >= 2)
            {
                fake_element.faces[0] = (face *) add_face(face(fake_element.nodes[0], fake_element.nodes[1], fake_element.nodes[2]), faces);
                fake_element.faces[1] = (face *) add_face(face(fake_element.nodes[0], fake_element.nodes[1], fake_element.nodes[3]), faces);
                fake_element.faces[2] = (face *) add_face(face(fake_element.nodes[0], fake_element.nodes[2], fake_element.nodes[3]), faces);
                fake_element.faces[3] = (face *) add_face(face(fake_element.nodes[1], fake_element.nodes[2], fake_element.nodes[3]), faces);
            }

            push_back_wrapper(fes, fake_element);
        }
        else if(type_of_elem == MSH_TRI_3)
        {
            map<phys_id, phys_area>::iterator ph = phys.find(phys_id(MSH_TRI_3, phys_num));
            if(ph != phys.end())
                fake_triangle.phys = &(ph->second);
            else
            {
                cout << "Error in " << __FILE__ << ":" << __LINE__
                     << " while reading file " << gmsh_filename << endl;
                cout << "Can`t detect physical id " << phys_num << " in 2 - MSH_TRI_3 (triangle)" << endl;
                return false;
            }

            size_t bound_type = fake_triangle.phys->type_of_bounds;
            if(bound_type != 1 && bound_type != 2)
            {
                cout << "Error: unaccounted bound, breaking..." << endl;
                return false;
            }

            for(size_t j = 0; j < 3; j++)
                gmsh_file >> local_nodes_tr[j];
            sort(local_nodes_tr.begin(), local_nodes_tr.end());
            for(size_t j = 0; j < 3; j++)
                fake_triangle.nodes[j] = & nodes[local_nodes_tr[j] - 1];

            fake_triangle.edges[0] = (edge *) add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[1]), edges);
            fake_triangle.edges[1] = (edge *) add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[2]), edges);
            fake_triangle.edges[2] = (edge *) add_edge(edge(fake_triangle.nodes[1], fake_triangle.nodes[2]), edges);
            if(config.basis.order >= 2)
                fake_triangle.faces = (face *) add_face(face(fake_triangle.nodes[0], fake_triangle.nodes[1], fake_triangle.nodes[2]), faces);
            if(config.boundary_enabled && bound_type == 1)
            {
                add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[1]), edges_surf_temp);
                add_edge(edge(fake_triangle.nodes[0], fake_triangle.nodes[2]), edges_surf_temp);
                add_edge(edge(fake_triangle.nodes[1], fake_triangle.nodes[2]), edges_surf_temp);
                if(config.basis.order >= 2)
                    add_face(face(fake_triangle.nodes[0], fake_triangle.nodes[1], fake_triangle.nodes[2]), faces_surf_temp);
            }
            if(config.boundary_enabled)
                push_back_wrapper(trs_full, triangle_full(fake_triangle));
            else
                push_back_wrapper(trs_base, fake_triangle);
        }
        else if(type_of_elem == MSH_LIN_2)
        {
            map<phys_id, phys_area>::iterator ph = phys.find(phys_id(MSH_LIN_2, phys_num));
            if(ph != phys.end())
                fake_edge_src.phys = &(ph->second);
            else
            {
                cout << "Error in " << __FILE__ << ":" << __LINE__
                     << " while reading file " << gmsh_filename << endl;
                cout << "Can`t detect physical id " << phys_num << " in 1 - MSH_LIN_2 (edge)" << endl;
                return false;
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
            push_back_wrapper(edges_src, fake_edge_src);
        }
        else
        {
            getline(gmsh_file, line);
        }
    }

    gmsh_file.close();

#if defined USE_CXX11
    fes.shrink_to_fit();
    trs_base.shrink_to_fit();
    trs_full.shrink_to_fit();
    pss.shrink_to_fit();
    edges_src.shrink_to_fit();
#else
    vector<finite_element>(fes).swap(fes);
    vector<triangle_base>(trs_base).swap(trs_base);
    vector<triangle_full>(trs_full).swap(trs_full);
    vector<pair<point, cvector3> >(pss).swap(pss);
    vector<edge_src>(edges_src).swap(edges_src);
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

    // Индексируем грани
    vector<face *> faces_ind;
    faces_ind.resize(faces.size());
    size_t index_face = 0;
    /// WARNING: const_cast для элемента set'а!
    for(set<face>::iterator i = faces.begin(); i != faces.end(); ++i)
    {
        assert(config.basis.order >= 2);
        face * fc = const_cast<face *>(&(*i));
        faces_ind[fc->num] = fc;
        fc->num = index_face++;
    }

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
        cout << "Error: 0 tetrahedrons detected, breaking..." << endl;
        return false;
    }
    for(size_t i = 0; i < fes.size(); i++)
    {
        show_progress("tetrahedrons", i, fes.size());
        // Переносим ребра и грани
        for(size_t j = 0; j < 6; j++)
            fes[i].edges[j] = edges_ind[(size_t)fes[i].edges[j]];
        if(config.basis.order >= 2)
            for(size_t j = 0; j < 4; j++)
                fes[i].faces[j] = faces_ind[(size_t)fes[i].faces[j]];
        // Инициализируем
        fes[i].init(& config.basis);
    }

    // Разбираемся с треугольниками
    trs.reserve(trs_base.size() + trs_full.size());
    for(size_t i = 0; i < trs_base.size(); i++)
        push_back_wrapper(trs, &trs_base[i]);
    for(size_t i = 0; i < trs_full.size(); i++)
        push_back_wrapper(trs, &trs_full[i]);
    if(trs.size() == 0)
    {
        cout << "Error: 0 triangles detected, breaking..." << endl;
        return false;
    }
    for(size_t i = 0; i < trs.size(); i++)
    {
        show_progress("triangles", i, trs.size());
        for(size_t j = 0; j < 3; j++)
            trs[i]->edges[j] = edges_ind[(size_t)trs[i]->edges[j]];
        if(config.basis.order >= 2)
            trs[i]->faces = faces_ind[(size_t)trs[i]->faces];

        // Заполняем степени свободы для первых краевых
        if(trs[i]->phys->type_of_bounds == 1)
        {
            array_t<size_t> dof(config.basis.tr_bf_num);
            array_t<size_t> dof_surf(config.basis.tr_bf_num);
            array_t<size_t> ker_dof(config.basis.tr_ker_bf_num);
            // Первый неполный
            for(size_t j = 0; j < 3; j++)
                dof[j] = trs[i]->edges[j]->num;
            for(size_t j = 0; j < 3; j++)
                ker_dof[j] = trs[i]->nodes[j]->num;
            if(config.boundary_enabled)
                for(size_t j = 0; j < 3; j++)
                    dof_surf[j] = edges_surf_temp.find(* trs[i]->edges[j])->num;
            // Первый полный
            if(config.basis.order >= 2 || config.basis.type == 2)
            {
                for(size_t j = 0; j < 3; j++)
                    dof[j + 3] = trs[i]->edges[j]->num + edges.size();
                for(size_t j = 0; j < 3; j++)
                    ker_dof[j + 3] = trs[i]->edges[j]->num + nodes.size();
                if(config.boundary_enabled)
                    for(size_t j = 0; j < 3; j++)
                        dof_surf[j + 3] = edges_surf_temp.find(* trs[i]->edges[j])->num + edges_surf_temp.size();
            }
            // Второй неполный
            if(config.basis.order >= 2)
            {
                dof[6] = trs[i]->faces->num + 2 * edges.size();
                dof[7] = trs[i]->faces->num + 2 * edges.size() + faces.size();
                if(config.boundary_enabled)
                {
                    dof_surf[6] = faces_surf_temp.find(* trs[i]->faces)->num + 2 * edges_surf_temp.size();
                    dof_surf[7] = faces_surf_temp.find(* trs[i]->faces)->num + 2 * edges_surf_temp.size() + faces_surf_temp.size();
                }
            }
            // Второй полный
            if(config.basis.order > 2 || (config.basis.type == 2 && config.basis.order == 2))
            {
                dof[8] = trs[i]->faces->num + 2 * edges.size() + 2 * faces.size();
                for(size_t j = 0; j < 3; j++)
                    dof[j + 9] = trs[i]->edges[j]->num + 2 * edges.size() + 3 * faces.size();
                ker_dof[6] = trs[i]->faces->num + edges.size() + nodes.size();
                for(size_t j = 0; j < 3; j++)
                    ker_dof[j + 7] = trs[i]->edges[j]->num + edges.size() + faces.size() + nodes.size();
                if(config.boundary_enabled)
                {
                    dof_surf[8] = faces_surf_temp.find(* trs[i]->faces)->num + 2 * edges_surf_temp.size() + 2 * faces_surf_temp.size();
                    for(size_t j = 0; j < 3; j++)
                        dof_surf[j + 9] = edges_surf_temp.find(* trs[i]->edges[j])->num + 2 * edges_surf_temp.size() + 3 * faces_surf_temp.size();
                }
            }
            // Теперь все собираем
            if(config.boundary_enabled)
                for(size_t j = 0; j < config.basis.tr_bf_num; j++)
                    global_to_local[dof[j]] = dof_surf[j];
            else
                for(size_t j = 0; j < config.basis.tr_bf_num; j++)
                    dof_first.insert(dof[j]);
            for(size_t j = 0; j < config.basis.tr_ker_bf_num; j++)
                ker_dof_first.insert(ker_dof[j]);
        }

        // Инициализируем
        trs[i]->init(& config.basis);
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
              min_coord[2], max_coord[2], fes);

#if defined VFEM_USE_PML
    input_pml();
#endif

    // Подсчитаем общее число степеней свободы
    dof_num = ker_dof_num = 0;
    if(config.basis.order == 1)
    {
        if(config.basis.type == 1)
        {
            dof_num = edges.size();
            ker_dof_num = nodes.size();
        }
        else if(config.basis.type == 2)
        {
            dof_num = 2 * edges.size();
            ker_dof_num = nodes.size() + edges.size();
        }
    }
    else if(config.basis.order == 2)
    {
        if(config.basis.type == 1)
        {
            dof_num = 2 * edges.size() + 2 * faces.size();
            ker_dof_num = nodes.size() + edges.size();
        }
        else if(config.basis.type == 2)
        {
            dof_num = 3 * edges.size() + 3 * faces.size();
            ker_dof_num = nodes.size() + 2 * edges.size() + faces.size();
        }
    }

    cout << "Statistics:" << endl;
    cout << " # Tetrehedrons:    " << fes.size() << endl;
    cout << " # Triangles:       " << trs.size() << endl;
    cout << " # Nodes:           " << nodes.size() << endl;
    cout << " # Edges:           " << edges.size() << endl;
    if(config.basis.order >= 2)
        cout << " # Faces:           " << faces.size() << endl;
    cout << " # Point sources:   " << pss.size() << endl;
    cout << " # Edge sources:    " << edges_src.size() << endl;
    cout << " # Main SLAE size:  " << dof_num << endl;
    cout << " # Ker. SLAE size:  " << ker_dof_num << endl;
    if(config.boundary_enabled)
        cout << " # Surf. SLAE size: " << global_to_local.size() << endl;

    return true;
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
            cvector3 s(0, 0, 0);
            if(flag[0] || flag[1] || flag[2])
                s = get_s(gauss_point_global, fefe, & phys_pml);

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
    phys_pml = config.phys_pml;

    if(phys_pml.params.size() == 0)
        return;

    cout << " > Building PML-bound coordinates ..." << endl;

    // Границы не PML точек
    double x0 = DBL_MAX;
    double x1 = -DBL_MAX;
    double y0 = DBL_MAX;
    double y1 = -DBL_MAX;
    double z0 = DBL_MAX;
    double z1 = -DBL_MAX;
    map<point *, pair<cpoint, finite_element *> > pml_nodes_cache;

    for(size_t i = 0; i < fes.size(); i++)
    {
        show_progress("scanned tetrahedrons", i, fes.size());
        finite_element * fes_i = &(fes[i]);
        if(is_pml(fes_i->barycenter, fes_i, & phys_pml))
        {
            for(size_t j = 0; j < 4; j++)
                pml_nodes_cache[fes_i->nodes[j]] = make_pair(cpoint(), fes_i);
        }
        else
        {
            for(size_t j = 0; j < 4; j++)
            {
                point * nodes_j = fes_i->nodes[j];
                if(nodes_j->x > x1) x1 = nodes_j->x;
                if(nodes_j->y > y1) y1 = nodes_j->y;
                if(nodes_j->z > z1) z1 = nodes_j->z;
                if(nodes_j->x < x0) x0 = nodes_j->x;
                if(nodes_j->y < y0) y0 = nodes_j->y;
                if(nodes_j->z < z0) z0 = nodes_j->z;
            }
        }
    }
    double big_num = DBL_MAX * 0.999;
    if(phys_pml.x0 >= big_num)  phys_pml.x0 = x0;
    if(phys_pml.x1 <= -big_num) phys_pml.x1 = x1;
    if(phys_pml.y0 >= big_num)  phys_pml.y0 = y0;
    if(phys_pml.y1 <= -big_num) phys_pml.y1 = y1;
    if(phys_pml.z0 >= big_num)  phys_pml.z0 = z0;
    if(phys_pml.z1 <= -big_num) phys_pml.z1 = z1;
    cout << "  non-PML dimension: (" << phys_pml.x0 << "," << phys_pml.x1 << ")x("
                                     << phys_pml.y0 << "," << phys_pml.y1 << ")x("
                                     << phys_pml.z0 << "," << phys_pml.z1 << ")" << endl;

    size_t i = 0;
    for(map<point *, pair<cpoint, finite_element *> >::iterator it = pml_nodes_cache.begin(); it != pml_nodes_cache.end(); ++it)
    {
        show_progress("replaced points", i++, pml_nodes_cache.size());

        point * p = it->first;
        pair<cpoint, finite_element *> * cp = &(it->second);
        finite_element * fefe = cp->second;
        cp->first = convert_point_to_pml(p, fefe);
    }

    for(size_t i = 0; i < fes.size(); i++)
    {
        show_progress("re-init tetrahedrons", i, fes.size());
        cpoint cp[4];
        size_t ph_curr = fes[i].phys->gmsh_num;
        if(is_pml(fes[i].barycenter, & (fes[i]), & phys_pml))
        {
            for(size_t j = 0; j < 4; j++)
            {
                map<point *, pair<cpoint, finite_element *> >::iterator it = pml_nodes_cache.find(fes[i].nodes[j]);
                // Если есть в кэше, то возьмем из него
                if(it->second.second->phys->gmsh_num == ph_curr)
                    cp[j] = it->second.first;
                // А иначе рассчитаем
                else
                    cp[j] = convert_point_to_pml(fes[i].nodes[j], &(fes[i]));
            }
            fes[i].init_pml(get_s, & phys_pml, cp);
        }
    }
}
#endif
