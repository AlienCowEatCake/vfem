#if !defined(EDGE_H_INCLUDED)
#define EDGE_H_INCLUDED

#include "../common/common.h"
#include "../vfem/phys.h"

// Структура ребро
struct edge : public edge_basic<point>
{
    edge(point * beg_ = NULL, point * end_ = NULL, size_t num_ = 0);
    edge(point & beg_, point & end_, size_t num_ = 0);
    double length();
};

// Структура ребро с источником
struct edge_src : public edge
{
    double direction;   // Направление
    phys_area * phys;   // Физическая область
    edge * edge_main;   // Основное ребро
    edge_src();
    edge_src(const edge_src & f);
    const phys_area & get_phys_area() const;
    edge_src & operator = (const edge_src & other);
};

#endif // EDGE_H_INCLUDED
