#include "edge.h"

// ============================================================================

edge::edge(point * beg_, point * end_, size_t num_)
{
    nodes[0] = beg_;
    nodes[1] = end_;
    num = num_;
}

edge::edge(point & beg_, point & end_, size_t num_)
{
    nodes[0] = & beg_;
    nodes[1] = & end_;
    num = num_;
}

double edge::length()
{
    double dx = nodes[1]->x - nodes[0]->x;
    double dy = nodes[1]->y - nodes[0]->y;
    double dz = nodes[1]->z - nodes[0]->z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

// ============================================================================

edge_src::edge_src()
{
    direction = 1.0;
    phys = NULL;
    edge_main = NULL;
}

edge_src::edge_src(const edge_src & f) : edge(f.nodes[0], f.nodes[1], f.num)
{
    direction = f.direction;
    phys = f.phys;
    edge_main = f.edge_main;
}

const phys_area & edge_src::get_phys_area() const
{
    assert(phys != NULL);
    return * phys;
}

edge_src & edge_src::operator = (const edge_src & other)
{
    if(this != & other)
    {
        this->nodes[0] = other.nodes[0];
        this->nodes[1] = other.nodes[1];
        this->num = other.num;
        this->direction = other.direction;
        this->phys = other.phys;
        this->edge_main = other.edge_main;
    }
    return * this;
}

// ============================================================================
