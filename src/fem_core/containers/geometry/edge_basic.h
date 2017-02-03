#if !defined(CONTAINERS_GEOMETRY_EDGE_BASIC_H_INCLUDED)
#define CONTAINERS_GEOMETRY_EDGE_BASIC_H_INCLUDED

#include <cassert>
#include <cstddef>

namespace fem_core { namespace containers { namespace geometry {

/**
 * @brief Структура ребро из двух узлов
 */
template<typename point>
struct edge_basic
{
    point * nodes[2];   ///< Узлы ребра
    std::size_t num;    ///< Номер ребра

    edge_basic(const edge_basic & f)
        : num(f.num)
    {
        nodes[0] = f.nodes[0];
        nodes[1] = f.nodes[1];
    }

    edge_basic(point * beg_ = NULL, point * end_ = NULL, std::size_t num_ = 0)
        : num(num_)
    {
        nodes[0] = beg_;
        nodes[1] = end_;
    }

    edge_basic(point & beg_, point & end_, std::size_t num_ = 0)
        : num(num_)
    {
        nodes[0] = & beg_;
        nodes[1] = & end_;
    }

    point & operator [] (std::size_t i)
    {
        assert(i < 2);
        return * nodes[i];
    }

    const point & operator [] (std::size_t i) const
    {
        assert(i < 2);
        return * nodes[i];
    }

    bool operator < (const edge_basic & t) const
    {
        assert(nodes[0] && nodes[1]);
        assert(t.nodes[0] && t.nodes[1]);
        return nodes[0]->num < t.nodes[0]->num || (!(t.nodes[0]->num < nodes[0]->num) && nodes[1]->num < t.nodes[1]->num);
    }

    bool operator == (const edge_basic & t) const
    {
        assert(nodes[0] && nodes[1]);
        assert(t.nodes[0] && t.nodes[1]);
        return nodes[0]->num == t.nodes[0]->num && nodes[1]->num == t.nodes[1]->num;
    }

    edge_basic & operator = (const edge_basic & other)
    {
        if(this != & other)
        {
            this->nodes[0] = other.nodes[0];
            this->nodes[1] = other.nodes[1];
            this->num = other.num;
        }
        return * this;
    }
};

}}} // namespace fem_core::containers::geometry

#endif // CONTAINERS_GEOMETRY_EDGE_BASIC_H_INCLUDED
