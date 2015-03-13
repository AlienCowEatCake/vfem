#include "face.h"

bool face::operator < (const face & t) const
{
    assert(nodes[0] && nodes[1] && nodes[2]);
    assert(t.nodes[0] && t.nodes[1] && t.nodes[2]);
    if(nodes[0] == t.nodes[0])
    {
        if(nodes[1] == t.nodes[1])
            return nodes[2] < t.nodes[2];
        return nodes[1] < t.nodes[1];
    }
    return nodes[0] < t.nodes[0];
}

bool face::operator == (const face & t) const
{
    assert(nodes[0] && nodes[1] && nodes[2]);
    assert(t.nodes[0] && t.nodes[1] && t.nodes[2]);
    return nodes[0]->num == t.nodes[0]->num &&
           nodes[1]->num == t.nodes[1]->num &&
           nodes[2]->num == t.nodes[2]->num;
}

face::face()
{
    nodes[0] = nodes[1] = nodes[2] = NULL;
    num = 0;
}

face::face(const face & f)
{
    for(size_t i = 0; i < 3; i++)
        nodes[i] = f.nodes[i];
    num = f.num;
}

face::face(point * p1, point * p2, point * p3)
{
    nodes[0] = p1;
    nodes[1] = p2;
    nodes[2] = p3;
    num = 0;
}

face::face(point & p1, point & p2, point & p3)
{
    nodes[0] = & p1;
    nodes[1] = & p2;
    nodes[2] = & p3;
    num = 0;
}

face::face(point * p1, point * p2, point * p3, size_t num_)
{
    nodes[0] = p1;
    nodes[1] = p2;
    nodes[2] = p3;
    num = num_;
}

face::face(point & p1, point & p2, point & p3, size_t num_)
{
    nodes[0] = & p1;
    nodes[1] = & p2;
    nodes[2] = & p3;
    num = num_;
}

point & face::operator [] (size_t i)
{
    assert(i < 3);
    return * nodes[i];
}

point face::operator [] (size_t i) const
{
    assert(i < 3);
    return * nodes[i];
}

face & face::operator = (const face & other)
{
    if(this != & other)
    {
        this->nodes[0] = other.nodes[0];
        this->nodes[1] = other.nodes[1];
        this->nodes[2] = other.nodes[2];
        this->num = other.num;
    }
    return * this;
}

ostream & operator << (ostream & os, const face & a)
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
    os << ", ";
    if(a.nodes[2])
        os << * a.nodes[2];
    else
        os << "NULL";
    os << " >";
    return os;
}
