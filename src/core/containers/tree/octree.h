#if !defined(CONTAINERS_TREE_OCTREE_H_INCLUDED)
#define CONTAINERS_TREE_OCTREE_H_INCLUDED

#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>

namespace core { namespace containers { namespace tree {

// *************************************************************************************************

template<class element_type>
class octree
{
public:
    void make(double x0, double x1, double y0, double y1, double z0, double z1, std::vector<element_type> & elements); //габариты области и список элементов
    element_type * find(double x, double y, double z) const;
    void clear();

private:

    class octree_node
    {
    public:
        octree_node();
        ~octree_node();
        void clear();
        void init_node(double x0, double x1, double y0, double y1, double z0, double z1, std::size_t split_constant, std::size_t level_barier);
        void add_element(std::vector<element_type *> & added_elements);
        bool add_element(element_type * added_element);
        element_type * find(double x, double y, double z) const;

        std::size_t level;

    private:
        bool point_in_node(double x, double y, double z) const;
        bool split_condition();
        void split();

        std::vector<element_type *> elements;
        octree_node * sub_nodes;
        double x0, x1, y0, y1, z0, z1;
        double eps_x, eps_y, eps_z;
        std::size_t split_constant;
        std::size_t level_barier;
    };

    octree_node root;
};

// *************************************************************************************************

template<class element_type>
octree<element_type>::octree_node::octree_node()
{
    sub_nodes = NULL;
}

template<class element_type>
octree<element_type>::octree_node::~octree_node()
{
    clear();
}

template<class element_type>
void octree<element_type>::octree_node::clear()
{
    if(sub_nodes)
    {
        for(std::size_t i = 0; i < 8; i++)
            sub_nodes[i].clear();
        delete [] sub_nodes;
        sub_nodes = NULL;
    }
    elements.clear();
}

template<class element_type>
void octree<element_type>::octree_node::init_node(double x0, double x1, double y0, double y1, double z0, double z1, std::size_t split_constant, std::size_t level_barier)
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

template<class element_type>
bool octree<element_type>::octree_node::point_in_node(double x, double y, double z) const
{
    if(x >= x0 - eps_x && x <= x1 + eps_x &&
       y >= y0 - eps_y && y <= y1 + eps_y &&
       z >= z0 - eps_z && z <= z1 + eps_z)
        return true;
    return false;
}

template<class element_type>
bool octree<element_type>::octree_node::split_condition()
{
    if(elements.size() + 1 == split_constant && fabs(x1 - x0) > eps_x && fabs(y1 - y0) > eps_y && fabs(z1 - z0) > eps_z && level <= level_barier)
        return true;
    return false;
}

template<class element_type>
bool octree<element_type>::octree_node::add_element(element_type * added_element)
{
    if(added_element->in_cube(x0, x1, y0, y1, z0, z1))
    {
        if(!sub_nodes)  //если узел терминальный (последний)
        {
            if(split_condition())   //если узел после добавляения перестанет быть терминальным
            {
                elements.push_back(added_element);
                split();
            }
            else   //если после добавления останется треминальным
                elements.push_back(added_element);
        }
        else   //если узел не терминальный (есть поддеревья)
        {
            for(std::size_t i = 0; i < 8; i++)
                sub_nodes[i].add_element(added_element);
        }
        return true;
    }
    return false;
}

template<class element_type>
void octree<element_type>::octree_node::add_element(std::vector<element_type *> & added_elements)
{
    for(std::size_t i = 0 ; i < added_elements.size(); i++)
        add_element(added_elements[i]);
}

template<class element_type>
void octree<element_type>::octree_node::split()
{
    sub_nodes = new octree_node [8];

    double dx = (x1 - x0) / 2.0, dy = (y1 - y0) / 2.0, dz = (z1 - z0) / 2.0;

    sub_nodes[0].init_node(x0, x0+dx, y0, y0+dy, z0, z0+dz, split_constant, level_barier);
    sub_nodes[1].init_node(x0+dx, x1, y0, y0+dy, z0, z0+dz, split_constant, level_barier);
    sub_nodes[2].init_node(x0, x0+dx, y0+dy, y1, z0, z0+dz, split_constant, level_barier);
    sub_nodes[3].init_node(x0+dx, x1, y0+dy, y1, z0, z0+dz, split_constant, level_barier);

    sub_nodes[4].init_node(x0, x0+dx, y0, y0+dy, z0+dz, z1, split_constant, level_barier);
    sub_nodes[5].init_node(x0+dx, x1, y0, y0+dy, z0+dz, z1, split_constant, level_barier);
    sub_nodes[6].init_node(x0, x0+dx, y0+dy, y1, z0+dz, z1, split_constant, level_barier);
    sub_nodes[7].init_node(x0+dx, x1, y0+dy, y1, z0+dz, z1, split_constant, level_barier);

    for(std::size_t i = 0; i < 8; i++)
        sub_nodes[i].level = level + 1;

    add_element(elements);
    elements.clear();
}

template<class element_type>
element_type * octree<element_type>::octree_node::find(double x, double y, double z) const
{
    if(point_in_node(x, y, z))
    {
        if(!sub_nodes)
        {
            std::size_t elements_num = elements.size();
            for(std::size_t i = 0; i < elements_num; i++)
                if(elements[i]->inside(x, y, z))
                    return elements[i];
        }
        else
        {
            for(std::size_t i = 0; i < 8; i++)
            {
                element_type * finded = sub_nodes[i].find(x, y, z);
                if(finded != NULL)
                    return finded;
            }
        }
    }
    return NULL;
}

template<class element_type>
void octree<element_type>::make(double x0, double x1, double y0, double y1, double z0, double z1, std::vector<element_type> & elements)
{
    std::size_t total_num = elements.size();
    std::size_t level_barier = static_cast<std::size_t>(std::log(static_cast<double>(total_num) + 1) / std::log(8.0)) + 1;
    std::size_t split_constant = static_cast<std::size_t>(std::pow(static_cast<double>(total_num), 1.0 / 8.0) * 1000.0);
    root.init_node(x0, x1, y0, y1, z0, z1, split_constant, level_barier);
    root.level = 0;
    std::vector<element_type *> el_v(total_num);
    for(std::size_t i = 0; i < total_num; i++)
        el_v[i] = &elements[i];
    root.add_element(el_v);
}

template<class element_type>
element_type * octree<element_type>::find(double x, double y, double z) const
{
    return root.find(x, y, z);
}

template<class element_type>
void octree<element_type>::clear()
{
    root.clear();
}

// *************************************************************************************************

}}} // core::containers::tree

#endif // CONTAINERS_TREE_OCTREE_H_INCLUDED
