#if !defined(CONTAINERS_TREE_QUADTREE_H_INCLUDED)
#define CONTAINERS_TREE_QUADTREE_H_INCLUDED

#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>

namespace core { namespace containers { namespace tree {

// *************************************************************************************************

template<class element_type>
class quadtree
{
public:
    void make(double x0, double x1, double y0, double y1, std::vector<element_type> & elements); //габариты области и список элементов
    element_type * find(double x, double y) const;
    void clear();

private:

    class quadtree_node
    {
    public:
        quadtree_node();
        ~quadtree_node();
        void clear();
        void init_node(double x0, double x1, double y0, double y1, size_t split_constant, size_t level_barier);
        void add_element(std::vector<element_type *> & added_elements);
        bool add_element(element_type * added_element);
        element_type * find(double x, double y) const;

        size_t level;

    private:
        bool point_in_node(double x, double y) const;
        bool split_condition();
        void split();

        std::vector<element_type *> elements;
        quadtree_node * sub_nodes;
        double x0, x1, y0, y1;
        double eps_x, eps_y;
        size_t split_constant;
        size_t level_barier;
    };

    quadtree_node root;
};

// ============================================================================

template<class element_type>
quadtree<element_type>::quadtree_node::quadtree_node()
{
    sub_nodes = NULL;
}

template<class element_type>
quadtree<element_type>::quadtree_node::~quadtree_node()
{
    clear();
}

template<class element_type>
void quadtree<element_type>::quadtree_node::clear()
{
    if(sub_nodes)
    {
        for(size_t i = 0; i < 4; i++)
            sub_nodes[i].clear();
        delete [] sub_nodes;
        sub_nodes = NULL;
    }
    elements.clear();
}

template<class element_type>
void quadtree<element_type>::quadtree_node::init_node(double x0, double x1, double y0, double y1, size_t split_constant, size_t level_barier)
{
    this->x0 = x0;
    this->x1 = x1;
    this->y0 = y0;
    this->y1 = y1;
    this->split_constant = split_constant;
    this->level_barier = level_barier;

    eps_x = 1e-12 * (x0 > x1 ? x0 : x1);
    eps_y = 1e-12 * (y0 > y1 ? y0 : y1);
}

template<class element_type>
bool quadtree<element_type>::quadtree_node::point_in_node(double x, double y) const
{
    if(x >= x0 - eps_x && x <= x1 + eps_x &&
       y >= y0 - eps_y && y <= y1 + eps_y)
        return true;
    return false;
}

template<class element_type>
bool quadtree<element_type>::quadtree_node::split_condition()
{
    if(elements.size() + 1 == split_constant && fabs(x1 - x0) > eps_x && fabs(y1 - y0) > eps_y && level <= level_barier)
        return true;
    return false;
}

template<class element_type>
bool quadtree<element_type>::quadtree_node::add_element(element_type * added_element)
{
    if(added_element->in_cube(x0, x1, y0, y1))
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
            for(size_t i = 0; i < 4; i++)
                sub_nodes[i].add_element(added_element);
        }
        return true;
    }
    return false;
}

template<class element_type>
void quadtree<element_type>::quadtree_node::add_element(std::vector<element_type *> & added_elements)
{
    for(size_t i = 0 ; i < added_elements.size(); i++)
        add_element(added_elements[i]);
}

template<class element_type>
void quadtree<element_type>::quadtree_node::split()
{
    sub_nodes = new quadtree_node [4];

    double dx = (x1 - x0) / 2.0, dy = (y1 - y0) / 2.0;

    sub_nodes[0].init_node(x0, x0+dx, y0, y0+dy, split_constant, level_barier);
    sub_nodes[1].init_node(x0+dx, x1, y0, y0+dy, split_constant, level_barier);
    sub_nodes[2].init_node(x0, x0+dx, y0+dy, y1, split_constant, level_barier);
    sub_nodes[3].init_node(x0+dx, x1, y0+dy, y1, split_constant, level_barier);

    for(size_t i = 0; i < 4; i++)
        sub_nodes[i].level = level + 1;

    add_element(elements);
    elements.clear();
}

template<class element_type>
element_type * quadtree<element_type>::quadtree_node::find(double x, double y) const
{
    if(point_in_node(x, y))
    {
        if(!sub_nodes)
        {
            for(size_t i = 0; i < elements.size(); i++)
                if(elements[i]->inside(x, y))
                    return elements[i];
        }
        else
        {
            for(size_t i = 0; i < 4; i++)
            {
                element_type * finded = sub_nodes[i].find(x, y);
                if(finded != NULL)
                    return finded;
            }
        }
    }
    return NULL;
}

template<class element_type>
void quadtree<element_type>::make(double x0, double x1, double y0, double y1, std::vector<element_type> & elements)
{
    size_t total_num = elements.size();
    size_t level_barier = (size_t)(log((double)total_num + 1) / log(4.0)) + 1;
    size_t split_constant = (size_t)(pow((double)total_num, 1.0 / 4.0) * 100.0);
    root.init_node(x0, x1, y0, y1, split_constant, level_barier);
    root.level = 0;
    std::vector<element_type *> el_v(total_num);
    for(size_t i = 0; i < total_num; i++)
        el_v[i] = &elements[i];
    root.add_element(el_v);
}

template<class element_type>
element_type * quadtree<element_type>::find(double x, double y) const
{
    return root.find(x, y);
}

template<class element_type>
void quadtree<element_type>::clear()
{
    root.clear();
}

// *************************************************************************************************

}}} // core::containers::tree

#endif // CONTAINERS_TREE_QUADTREE_H_INCLUDED
