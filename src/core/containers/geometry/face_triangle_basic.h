#if !defined(CONTAINERS_GEOMETRY_FACE_TRIANGLE_BASIC_H_INCLUDED)
#define CONTAINERS_GEOMETRY_FACE_TRIANGLE_BASIC_H_INCLUDED

#include <cassert>
#include <cstddef>

namespace core { namespace containers { namespace geometry {

/**
 * @brief Структура треугольная грань из трех узлов
 */
template<typename point>
struct face_triangle_basic
{
    point * nodes[3];   ///< Узлы грани
    std::size_t num;    ///< Номер грани

    face_triangle_basic()
    {
        nodes[0] = nodes[1] = nodes[2] = NULL;
        num = 0;
    }

    face_triangle_basic(const face_triangle_basic & f)
    {
        for(size_t i = 0; i < 3; i++)
            nodes[i] = f.nodes[i];
        num = f.num;
    }

    face_triangle_basic(point * p1, point * p2, point * p3, std::size_t num_ = 0)
    {
        nodes[0] = p1;
        nodes[1] = p2;
        nodes[2] = p3;
        num = num_;
    }

    face_triangle_basic(point & p1, point & p2, point & p3, std::size_t num_ = 0)
    {
        nodes[0] = & p1;
        nodes[1] = & p2;
        nodes[2] = & p3;
        num = num_;
    }

    point & operator [] (size_t i)
    {
        assert(i < 3);
        return * nodes[i];
    }

    const point & operator [] (size_t i) const
    {
        assert(i < 3);
        return * nodes[i];
    }

    bool operator < (const face_triangle_basic & t) const
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

    bool operator == (const face_triangle_basic & t) const
    {
        assert(nodes[0] && nodes[1] && nodes[2]);
        assert(t.nodes[0] && t.nodes[1] && t.nodes[2]);
        return nodes[0]->num == t.nodes[0]->num &&
               nodes[1]->num == t.nodes[1]->num &&
               nodes[2]->num == t.nodes[2]->num;
    }

    face_triangle_basic & operator = (const face_triangle_basic & other)
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
};

}}} // namespace core::containers::geometry

#endif // CONTAINERS_GEOMETRY_FACE_TRIANGLE_BASIC_H_INCLUDED
