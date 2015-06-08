#include "edge.h"

bool edge::operator < (const edge & t) const
{
    assert(nodes[0] && nodes[1]);
    assert(t.nodes[0] && t.nodes[1]);
    return nodes[0]->num < t.nodes[0]->num || (!(t.nodes[0]->num < nodes[0]->num) && nodes[1]->num < t.nodes[1]->num);
}

bool edge::operator == (const edge & t) const
{
    assert(nodes[0] && nodes[1]);
    assert(t.nodes[0] && t.nodes[1]);
    return nodes[0]->num == t.nodes[0]->num && nodes[1]->num == t.nodes[1]->num;
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
    assert(i < 2);
    return * nodes[i];
}

const point & edge::operator [] (size_t i) const
{
    assert(i < 2);
    return * nodes[i];
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
