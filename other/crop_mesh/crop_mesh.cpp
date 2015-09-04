#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <cstdlib>

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
    gmsh_file << "$Elements\n" << eds.size() + trs.size() + tets.size() << '\n';
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

int main(int argc, char * argv[])
{
    if(argc < 4)
    {
#if defined _WIN32
        char * a = strrchr(argv[0], '\\');
#else
        char * a = strrchr(argv[0], '/');
#endif
        char * progname = a ? a + 1 : argv[0];
        cerr << "Usage: " << progname << " input_mesh output_mesh phys_tag [phys_tag ...]" << endl;
#if defined _WIN32
        system("pause");
#endif
        return 1;
    }

    cout << "Reading mesh ..." << endl;
    if(!read_gmsh(argv[1]))
    {
        cerr << "Error reading file " << argv[1] << endl;
#if defined _WIN32
        system("pause");
#endif
        return 1;
    }

    set<size_t> exclude_tags;
    for(int i = 3; i < argc; i++)
        exclude_tags.insert(atoi(argv[i]));

    cout << "Cropping mesh ..." << endl;
    set<size_t> used_nodes;

    cout << " * tetrahedrons ..." << endl;
    for(size_t i = 0; i < tets.size(); i++)
        if(exclude_tags.find(tets[i].tags[0]) != exclude_tags.end())
            tets.erase(tets.begin() + i--);
        else
            for(size_t j = 0; j < 4; j++)
                used_nodes.insert(tets[i].nodes[j]);

    cout << " * triangles ..." << endl;
    for(size_t i = 0; i < trs.size(); i++)
        if(exclude_tags.find(trs[i].tags[0]) != exclude_tags.end())
            trs.erase(trs.begin() + i--);
        else
            for(size_t j = 0; j < 3; j++)
                used_nodes.insert(trs[i].nodes[j]);

    cout << " * edges ..." << endl;
    for(size_t i = 0; i < eds.size(); i++)
        if(exclude_tags.find(eds[i].tags[0]) != exclude_tags.end())
            eds.erase(eds.begin() + i--);
        else
            for(size_t j = 0; j < 2; j++)
                used_nodes.insert(eds[i].nodes[j]);

    cout << " * nodes ..." << endl;
    for(size_t i = 0; i < nodes.size(); i++)
        if(used_nodes.find(nodes[i].num) == used_nodes.end())
            nodes.erase(nodes.begin() + i--);

    map<size_t, size_t> old_to_new;
    for(size_t i = 0; i < nodes.size(); i++)
        old_to_new[nodes[i].num] = nodes[i].num = i + 1;

    for(size_t i = 0; i < tets.size(); i++)
        for(size_t j = 0; j < 4; j++)
            tets[i].nodes[j] = old_to_new[tets[i].nodes[j]];
    for(size_t i = 0; i < trs.size(); i++)
        for(size_t j = 0; j < 3; j++)
            trs[i].nodes[j] = old_to_new[trs[i].nodes[j]];
    for(size_t i = 0; i < eds.size(); i++)
        for(size_t j = 0; j < 2; j++)
            eds[i].nodes[j] = old_to_new[eds[i].nodes[j]];

    cout << "Writing mesh ..." << endl;
    if(!write_gmsh(argv[2]))
    {
        cerr << "Error writing file " << argv[2] << endl;
#if defined _WIN32
        system("pause");
#endif
        return 1;
    }

    cout << "Done." << endl;
#if defined _WIN32
    system("pause");
#endif
    return 0;
}

