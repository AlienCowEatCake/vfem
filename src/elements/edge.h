#if !defined EDGE_H_INCLUDED
#define EDGE_H_INCLUDED

#include "../common/common.h"
#include "../geometry/point.h"
#include "../vfem/phys.h"

// Класс ребро
class edge
{
public:
    point * nodes[2];   // Узлы ребра
    size_t num;         // Номер ребра
    edge();
    edge(const edge & f);
    edge(point * beg_, point * end_);
    edge(point & beg_, point & end_);
    edge(point * beg_, point * end_, size_t num_);
    edge(point & beg_, point & end_, size_t num_);
    point & operator [] (size_t i);
    const point & operator [] (size_t i) const;
    bool operator < (const edge & t) const;
    bool operator == (const edge & t) const;
    edge & operator = (const edge & other);
    friend ostream & operator << (ostream & os, const edge & a);
};

#endif // EDGE_H_INCLUDED
