#if !defined OCTAL_TREE_H_INCLUDED
#define OCTAL_TREE_H_INCLUDED

#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>

using namespace std;

template<class element_type> class octal_tree_node
{
public:
    octal_tree_node();
    void init_node(double x0, double x1, double y0, double y1, double z0, double z1, size_t split_constant, size_t level_barier);
    void add_element(vector<element_type *> & added_elements);
    bool add_element(element_type * added_element);
    element_type * find(double x, double y, double z) const;

    size_t level;

private:
    bool point_in_node(double x, double y, double z) const;
    bool split_condition();
    void split();

    vector<element_type *> elements;
    bool node_is_terminal;
    octal_tree_node<element_type> * sub_nodes;
    double x0, x1, y0, y1, z0, z1;
    double eps_x, eps_y, eps_z;
    size_t elements_num;
    size_t split_constant;
    size_t level_barier;
};

template<class element_type> class octal_tree
{
public:
    void make(double x0, double x1, double y0, double y1, double z0, double z1, element_type * elements, size_t elements_num); //габариты области и список элементов
    element_type * find(double x, double y, double z) const;

private:
    octal_tree_node<element_type> root;
    size_t total_num;
    size_t split_constant;
};

// ============================================================================

template<class element_type> octal_tree_node<element_type>::octal_tree_node()
{
    node_is_terminal = true;
    elements_num = 0;
}

template<class element_type> void octal_tree_node<element_type>::init_node(double x0, double x1, double y0, double y1, double z0, double z1, size_t split_constant, size_t level_barier)
{
    this->x0 = x0;
    this->x1 = x1;
    this->y0 = y0;
    this->y1 = y1;
    this->z0 = z0;
    this->z1 = z1;
    this->split_constant = split_constant;
    this->level_barier = level_barier;

    eps_x = 1e-12 * (x0 > x1 ? x0 : x1);
    eps_y = 1e-12 * (y0 > y1 ? y0 : y1);
    eps_z = 1e-12 * (z0 > z1 ? z0 : z1);
}

template<class element_type> bool octal_tree_node<element_type>::point_in_node(double x, double y, double z) const
{
    if(x >= x0 - eps_x && x <= x1 + eps_x &&
       y >= y0 - eps_y && y <= y1 + eps_y &&
       z >= z0 - eps_z && z <= z1 + eps_z)
        return true;
    return false;
}

template<class element_type> bool octal_tree_node<element_type>::split_condition()
{
    if(elements_num + 1 == split_constant && fabs(x1 - x0) > eps_x && fabs(y1 - y0) > eps_y && fabs(z1 - z0) > eps_z && level <= level_barier)
        return true;
    return false;
}

template<class element_type> bool octal_tree_node<element_type>::add_element(element_type * added_element)
{
    if(added_element->inside_tree(x0, x1, y0, y1, z0, z1))
    {
        if(node_is_terminal)   //если узел терминальный (последний)
        {
            if(split_condition())   //если узел после добавляения перестанет быть терминальным
            {
                elements.push_back(added_element);
                split();
            }
            else   //если после добавления останется треминальным
                elements.push_back(added_element);
            elements_num++;
        }
        else   //если узел не терминальный (есть поддеревья)
        {
            elements_num++;
            for(size_t i = 0; i < 8; i++)
                sub_nodes[i].add_element(added_element);
        }
        return true;
    }
    return false;
}

template<class element_type> void octal_tree_node<element_type>::add_element(vector<element_type *> & added_elements)
{
    for(size_t i = 0 ; i < added_elements.size(); i++)
        add_element(added_elements[i]);
}

template<class element_type> void octal_tree_node<element_type>::split()
{
    node_is_terminal = false;
    sub_nodes = new octal_tree_node<element_type> [8];

    double dx = (x1 - x0) / 2.0, dy = (y1 - y0) / 2.0, dz = (z1 - z0) / 2.0;

    sub_nodes[0].init_node(x0, x0+dx, y0, y0+dy, z0, z0+dz, split_constant, level_barier);
    sub_nodes[1].init_node(x0+dx, x1, y0, y0+dy, z0, z0+dz, split_constant, level_barier);
    sub_nodes[2].init_node(x0, x0+dx, y0+dy, y1, z0, z0+dz, split_constant, level_barier);
    sub_nodes[3].init_node(x0+dx, x1, y0+dy, y1, z0, z0+dz, split_constant, level_barier);

    sub_nodes[4].init_node(x0, x0+dx, y0, y0+dy, z0+dz, z1, split_constant, level_barier);
    sub_nodes[5].init_node(x0+dx, x1, y0, y0+dy, z0+dz, z1, split_constant, level_barier);
    sub_nodes[6].init_node(x0, x0+dx, y0+dy, y1, z0+dz, z1, split_constant, level_barier);
    sub_nodes[7].init_node(x0+dx, x1, y0+dy, y1, z0+dz, z1, split_constant, level_barier);

    for(size_t i = 0; i < 8; i++)
        sub_nodes[i].level = level + 1;

    elements_num = 0;
    add_element(elements);
    elements.clear();
}

template<class element_type> element_type * octal_tree_node<element_type>::find(double x, double y, double z) const
{
    if(point_in_node(x, y, z))
    {
        if(node_is_terminal)
        {
            for(size_t i = 0; i < elements_num; i++)
                if (elements[i]->inside(x, y, z))
                    return elements[i];
        }
        else
        {
            for(size_t i = 0; i < 8; i++)
            {
                element_type * finded = sub_nodes[i].find(x, y, z);
                if(finded != NULL)
                    return finded;
            }
        }
    }
    return NULL;
}

template<class element_type> void octal_tree<element_type>::make(double x0, double x1, double y0, double y1, double z0, double z1, element_type * elements, size_t elements_num)
{
    total_num = elements_num;
    split_constant = (size_t)sqrt((double)total_num);
    size_t level_barier = (size_t)(log((double)total_num) / log(8.0)) + 5; //скорость поиска в дереве O(log(n)), где n - число элементов,
    root.init_node(x0, x1, y0, y1, z0, z1, split_constant, level_barier);
    root.level = 0;
    vector<element_type *> el_v;
    for(size_t i = 0; i < total_num; i++)
        el_v.push_back(&elements[i]);
    root.add_element(el_v);
}

template<class element_type> element_type * octal_tree<element_type>::find(double x, double y, double z) const
{
    return root.find(x, y, z);
}

#endif // OCTAL_TREE_H_INCLUDED
