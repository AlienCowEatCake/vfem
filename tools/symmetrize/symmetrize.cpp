#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
using namespace std;

struct node
{
    double coord[3];
    size_t num;
};

struct edge
{
    size_t nodes[2];
    size_t num;
    size_t tags_num;
    size_t tags[2];
};

struct triangle
{
    size_t nodes[3];
    size_t num;
    size_t tags_num;
    size_t tags[2];
};

struct tetrahedron
{
    size_t nodes[4];
    size_t num;
    size_t tags_num;
    size_t tags[2];
};

vector<edge> eds;
vector<node> nodes;
vector<triangle> trs;
vector<tetrahedron> tets;
map<size_t, size_t> new_nodes;

bool read_gmsh(string input)
{
    ifstream gmsh_file;
    gmsh_file.open(input.c_str(), ios::in);
    string line;

    do
        getline(gmsh_file, line);
    while(line.find("$Nodes") == string::npos && gmsh_file.good());

    if(!gmsh_file.good()) return false;

    size_t nodes_num;
    gmsh_file >> nodes_num;
    nodes.resize(nodes_num);

    for(size_t i = 0; i < nodes_num; i++)
    {
        gmsh_file >> nodes[i].num;
        for(size_t j = 0; j < 3; j++)
            gmsh_file >> nodes[i].coord[j];
    }

    do
        getline(gmsh_file, line);
    while(line.find("$Elements") == string::npos && gmsh_file.good());

    if(!gmsh_file.good()) return false;

    edge fake_ed;
    triangle fake_tr;
    tetrahedron fake_tet;

    size_t fes_num;
    gmsh_file >> fes_num;
    for(size_t i = 0; i < fes_num; i++)
    {
        size_t type_of_elem, fe_number;
        gmsh_file >> fe_number >> type_of_elem;
        if(type_of_elem == 1)
        {
            fake_ed.num = fe_number;
            gmsh_file >> fake_ed.tags_num;
            if(fake_ed.tags_num == 2)
                gmsh_file >> fake_ed.tags[0] >> fake_ed.tags[1];
            for(size_t j = 0; j < 2; j++)
                gmsh_file >> fake_ed.nodes[j];
            eds.push_back(fake_ed);
        }
        if(type_of_elem == 2)
        {
            fake_tr.num = fe_number;
            gmsh_file >> fake_tr.tags_num;
            if(fake_tr.tags_num == 2)
                gmsh_file >> fake_tr.tags[0] >> fake_tr.tags[1];
            for(size_t j = 0; j < 3; j++)
                gmsh_file >> fake_tr.nodes[j];
            trs.push_back(fake_tr);
        }
        else if(type_of_elem == 4)
        {
            fake_tet.num = fe_number;
            gmsh_file >> fake_tet.tags_num;
            if(fake_tet.tags_num == 2)
                gmsh_file >> fake_tet.tags[0] >> fake_tet.tags[1];
            for(size_t j = 0; j < 4; j++)
                gmsh_file >> fake_tet.nodes[j];
            tets.push_back(fake_tet);
        }
    }

    if(!gmsh_file.good()) return false;
    gmsh_file.close();
    return true;
}

bool write_gmsh(string output)
{
    ofstream gmsh_file;
    gmsh_file.open(output.c_str(), ios::out);
    if(!gmsh_file.good()) return false;

    //gmsh_file << scientific;
    //gmsh_file.precision(17);

    gmsh_file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";

    gmsh_file << "$Nodes\n" << nodes.size() << '\n';
    for(size_t i = 0; i < nodes.size(); i++)
    {
        gmsh_file << nodes[i].num << ' ';
        for(size_t j = 0; j < 3; j++)
            gmsh_file << nodes[i].coord[j] << ' ';
        gmsh_file << '\n';
    }
    gmsh_file << "$EndNodes\n";

    size_t curr_fe = 1;
    gmsh_file << "$Elements\n" << trs.size() + tets.size() << '\n';
    for(size_t i = 0; i < eds.size(); i++)
    {
        gmsh_file << curr_fe << " 1 " << eds[i].tags_num << ' ';
        for(size_t j = 0; j < 2; j++)
            gmsh_file << eds[i].tags[j] << ' ';
        for(size_t j = 0; j < 2; j++)
            gmsh_file << eds[i].nodes[j] << ' ';
        gmsh_file << '\n';
        curr_fe++;
    }
    for(size_t i = 0; i < trs.size(); i++)
    {
        gmsh_file << curr_fe << " 2 " << trs[i].tags_num << ' ';
        for(size_t j = 0; j < 2; j++)
            gmsh_file << trs[i].tags[j] << ' ';
        for(size_t j = 0; j < 3; j++)
            gmsh_file << trs[i].nodes[j] << ' ';
        gmsh_file << '\n';
        curr_fe++;
    }
    for(size_t i = 0; i < tets.size(); i++)
    {
        gmsh_file << curr_fe << " 4 " << tets[i].tags_num << ' ';
        for(size_t j = 0; j < 2; j++)
            gmsh_file << tets[i].tags[j] << ' ';
        for(size_t j = 0; j < 4; j++)
            gmsh_file << tets[i].nodes[j] << ' ';
        gmsh_file << '\n';
        curr_fe++;
    }
    gmsh_file << "$EndElements\n\n";

    if(!gmsh_file.good()) return false;
    gmsh_file.flush();
    gmsh_file.close();
    return true;
}

void symmetrize(string input, string output)
{
    if(!read_gmsh(input))
    {
        cerr << "Error reading file " << input << endl;
        return;
    }

    size_t nodes_size = nodes.size();
    for(size_t i = 0; i < nodes_size; i++)
    {
        node tmp_node = nodes[i];
        tmp_node.coord[0] *= -1.0;
        if(fabs(tmp_node.coord[0] - nodes[i].coord[0]) > 1.0e-15)
        {
            tmp_node.num = nodes.size() + 1;
            nodes.push_back(tmp_node);
        }
        new_nodes[i+1] = tmp_node.num;
    }

    size_t eds_size = eds.size();
    for(size_t i = 0; i < eds_size; i++)
    {
        edge tmp_edge = eds[i];
        // Warning: edges is rotated!
        tmp_edge.nodes[0] = new_nodes[tmp_edge.nodes[1]];
        tmp_edge.nodes[1] = new_nodes[tmp_edge.nodes[0]];
        eds.push_back(tmp_edge);
    }

    size_t trs_size = trs.size();
    for(size_t i = 0; i < trs_size; i++)
    {
        triangle tmp_triangle = trs[i];
        for(size_t j = 0; j < 3; j++)
            tmp_triangle.nodes[j] = new_nodes[tmp_triangle.nodes[j]];
        trs.push_back(tmp_triangle);
    }

    size_t tet_size = tets.size();
    for(size_t i = 0; i < tet_size; i++)
    {
        tetrahedron tmp_tetrahedron = tets[i];
        for(size_t j = 0; j < 4; j++)
            tmp_tetrahedron.nodes[j] = new_nodes[tmp_tetrahedron.nodes[j]];
        tets.push_back(tmp_tetrahedron);
    }

    if(!write_gmsh(output))
    {
        cerr << "Error writing file " << output << endl;
        return;
    }

    cout << "Done." << endl;
}

int main()
{
    symmetrize("area.msh", "area_symm.msh");
#if defined _WIN32
    system("pause");
#endif
    return 0;
}
