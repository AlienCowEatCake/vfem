#if !defined FACE_H_INCLUDED
#define FACE_H_INCLUDED

#include "../common/common.h"
#include "../geometry/point.h"

// Класс грань
class face
{
public:
    point * nodes[3];   // Узлы грани
    size_t num;         // Номер грани
    face();
    face(const face & f);
    face(point * p1, point * p2, point * p3);
    face(point & p1, point & p2, point & p3);
    face(point * p1, point * p2, point * p3, size_t num_);
    face(point & p1, point & p2, point & p3, size_t num_);
    point & operator [] (size_t i);
    const point & operator [] (size_t i) const;
    bool operator < (const face & t) const;
    bool operator == (const face & t) const;
    face & operator = (const face & other);
    friend ostream & operator << (ostream & os, const face & a);
};

#endif // FACE_H_INCLUDED
