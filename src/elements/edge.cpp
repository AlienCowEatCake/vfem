#include "edge.h"

// ============================================================================

bool edge::operator < (const edge & t) const
{
    if(nodes[0] && nodes[1])
        return nodes[0]->num < t.nodes[0]->num || (!(t.nodes[0]->num < nodes[0]->num) && nodes[1]->num < t.nodes[1]->num);
    else
        return false;
}

bool edge::operator == (const edge & t) const
{
    if(nodes[0] && nodes[1])
        return nodes[0]->num == t.nodes[0]->num && nodes[1]->num == t.nodes[1]->num;
    else
        return false;
}

edge::edge()
{
    nodes[0] = nodes[1] = NULL;
    num = 0;
}

edge::edge(const edge & f)
{
    nodes[0] = f.nodes[0];
    nodes[1] = f.nodes[1];
    num = f.num;
}

edge::edge(point * beg_, point * end_)
{
    nodes[0] = beg_;
    nodes[1] = end_;
    num = 0;
}

edge::edge(point & beg_, point & end_)
{
    nodes[0] = & beg_;
    nodes[1] = & end_;
    num = 0;
}

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

point & edge::operator [] (size_t i)
{
    switch(i)
    {
    case 0:
        return * nodes[0];
    case 1:
        return * nodes[1];
    default:
        cerr << "Error: Unknown index " << i << " at edge " << * this << endl;
        throw ADDRESSING_ERROR;
    }
}

point edge::operator [] (size_t i) const
{
    switch(i)
    {
    case 0:
        return * nodes[0];
    case 1:
        return * nodes[1];
    default:
        cerr << "Error: Unknown index " << i << " at edge " << * this << endl;
        throw ADDRESSING_ERROR;
    }
}

edge & edge::operator = (const edge & other)
{
    if(this != & other)
    {
        this->nodes[0] = other.nodes[0];
        this->nodes[1] = other.nodes[1];
        this->num = other.num;
    }
    return * this;
}

ostream & operator << (ostream & os, const edge & a)
{
    os << "< ";
    if(a.nodes[0])
        os << * a.nodes[0];
    else
        os << "NULL";
    os << ", ";
    if(a.nodes[1])
        os << * a.nodes[1];
    else
        os << "NULL";
    os << " >";
    return os;
}

// ============================================================================

edge_src::edge_src()
{
    direction = 1.0;
    phys = NULL;
    edge_main = NULL;
}

const phys_area & edge_src::get_phys_area() const
{
    if(!phys)
    {
        cerr << "Error: Null pointer at get_phys_area()" << endl;
        throw NULL_PTR_ERROR;
    }
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
