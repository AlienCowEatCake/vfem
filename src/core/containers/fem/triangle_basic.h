#if !defined(CONTAINERS_FEM_TRIANGLE_BASIC_H_INCLUDED)
#define CONTAINERS_FEM_TRIANGLE_BASIC_H_INCLUDED

#include <cassert>
#include <cstddef>

namespace core { namespace containers { namespace fem {

/**
 * @brief Класс треугольник
 */
template<typename point, typename edge, typename face, typename phys_area>
class triangle_basic
{
public:
    triangle_basic()
    {
        for(std::size_t i = 0; i < 3; i++)
            m_nodes[i] = NULL;
        for(std::size_t i = 0; i < 3; i++)
            m_edges[i] = NULL;
        m_faces = NULL;
        m_phys = NULL;
    }

    inline const point & get_node(std::size_t i) const
    {
        assert(i < 3);
        assert(m_nodes[i] != NULL);
        return (* m_nodes[i]);
    }

    inline point & get_node(std::size_t i)
    {
        assert(i < 3);
        assert(m_nodes[i] != NULL);
        return (* m_nodes[i]);
    }

    inline const point * get_node_ptr(std::size_t i) const
    {
        assert(i < 3);
        return m_nodes[i];
    }

    inline point * get_node_ptr(std::size_t i)
    {
        assert(i < 3);
        return m_nodes[i];
    }

    inline void set_node(std::size_t i, point * new_node)
    {
        assert(i < 3);
        m_nodes[i] = new_node;
    }

    inline void set_node(std::size_t i, point & new_node)
    {
        assert(i < 3);
        m_nodes[i] = & new_node;
    }

    inline const edge & get_edge(std::size_t i) const
    {
        assert(i < 3);
        assert(m_edges[i] != NULL);
        return (* m_edges[i]);
    }

    inline edge & get_edge(std::size_t i)
    {
        assert(i < 3);
        assert(m_edges[i] != NULL);
        return (* m_edges[i]);
    }

    inline const edge * get_edge_ptr(std::size_t i) const
    {
        assert(i < 3);
        return m_edges[i];
    }

    inline edge * get_edge_ptr(std::size_t i)
    {
        assert(i < 3);
        return m_edges[i];
    }

    inline void set_edge(std::size_t i, edge * new_edge)
    {
        assert(i < 3);
        m_edges[i] = new_edge;
    }

    inline void set_edge(std::size_t i, edge & new_edge)
    {
        assert(i < 3);
        m_edges[i] = & new_edge;
    }

    inline const face & get_face() const
    {
        assert(m_faces != NULL);
        return * m_faces;
    }

    inline face & get_face()
    {
        assert(m_faces != NULL);
        return * m_faces;
    }

    inline const face * get_face_ptr() const
    {
        return m_faces;
    }

    inline face * get_face_ptr()
    {
        return m_faces;
    }

    inline void set_face(face * new_face)
    {
        m_faces = new_face;
    }

    inline void set_face(face & new_face)
    {
        m_faces = & new_face;
    }

    inline const phys_area & get_phys_area() const
    {
        assert(m_phys != NULL);
        return * m_phys;
    }

    inline phys_area & get_phys_area()
    {
        assert(m_phys != NULL);
        return * m_phys;
    }

    inline const phys_area * get_phys_area_ptr() const
    {
        return m_phys;
    }

    inline phys_area * get_phys_area_ptr()
    {
        return m_phys;
    }

    inline void set_phys_area(phys_area * new_phys_area)
    {
        m_phys = new_phys_area;
    }

    inline void set_phys_area(phys_area & new_phys_area)
    {
        m_phys = & new_phys_area;
    }

private:

    point * m_nodes[3];   // Узлы
    edge * m_edges[3];    // Ребра
    face * m_faces;       // Грани
    phys_area * m_phys;   // Физическая область
};

}}} // namespace core::containers::fem

#endif // CONTAINERS_FEM_TRIANGLE_BASIC_H_INCLUDED
