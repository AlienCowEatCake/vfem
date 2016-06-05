#if !defined(CONTAINERS_FEM_TETRAHEDRON_BASIC_H_INCLUDED)
#define CONTAINERS_FEM_TETRAHEDRON_BASIC_H_INCLUDED

#include <cassert>
#include <cstddef>

namespace core { namespace containers { namespace fem {

/**
 * @brief Класс тетраэдр
 */
template<typename point, typename edge, typename face, typename phys_area>
class tetrahedron_basic
{
public:
    tetrahedron_basic()
    {
        for(std::size_t i = 0; i < 4; i++)
            m_nodes[i] = NULL;
        for(std::size_t i = 0; i < 6; i++)
            m_edges[i] = NULL;
        for(std::size_t i = 0; i < 4; i++)
            m_faces[i] = NULL;
        m_phys = NULL;
    }

    inline const point & get_node(std::size_t i) const
    {
        assert(i < 4);
        assert(m_nodes[i] != NULL);
        return (* m_nodes[i]);
    }

    inline point & get_node(std::size_t i)
    {
        assert(i < 4);
        assert(m_nodes[i] != NULL);
        return (* m_nodes[i]);
    }

    inline const point * get_node_ptr(std::size_t i) const
    {
        assert(i < 4);
        return m_nodes[i];
    }

    inline point * get_node_ptr(std::size_t i)
    {
        assert(i < 4);
        return m_nodes[i];
    }

    inline void set_node(std::size_t i, point * new_node)
    {
        assert(i < 4);
        m_nodes[i] = new_node;
    }

    inline void set_node(std::size_t i, point & new_node)
    {
        assert(i < 4);
        m_nodes[i] = & new_node;
    }

    inline const edge & get_edge(std::size_t i) const
    {
        assert(i < 6);
        assert(m_edges[i] != NULL);
        return (* m_edges[i]);
    }

    inline edge & get_edge(std::size_t i)
    {
        assert(i < 6);
        assert(m_edges[i] != NULL);
        return (* m_edges[i]);
    }

    inline const edge * get_edge_ptr(std::size_t i) const
    {
        assert(i < 6);
        return m_edges[i];
    }

    inline edge * get_edge_ptr(std::size_t i)
    {
        assert(i < 6);
        return m_edges[i];
    }

    inline void set_edge(std::size_t i, edge * new_edge)
    {
        assert(i < 6);
        m_edges[i] = new_edge;
    }

    inline void set_edge(std::size_t i, edge & new_edge)
    {
        assert(i < 6);
        m_edges[i] = & new_edge;
    }

    inline const face & get_face(std::size_t i) const
    {
        assert(i < 4);
        assert(m_faces != NULL);
        return (* m_faces[i]);
    }

    inline face & get_face(std::size_t i)
    {
        assert(i < 4);
        assert(m_faces != NULL);
        return (* m_faces[i]);
    }

    inline const face * get_face_ptr(std::size_t i) const
    {
        assert(i < 4);
        return m_faces[i];
    }

    inline face * get_face_ptr(std::size_t i)
    {
        assert(i < 4);
        return m_faces[i];
    }

    inline void set_face(std::size_t i, face * new_face)
    {
        assert(i < 4);
        m_faces[i] = new_face;
    }

    inline void set_face(std::size_t i, face & new_face)
    {
        assert(i < 4);
        m_faces[i] = & new_face;
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

    point * m_nodes[4];   // Узлы
    edge * m_edges[6];    // Ребра
    face * m_faces[4];    // Грани
    phys_area * m_phys;   // Физическая область
};

}}} // namespace core::containers::fem

#endif // CONTAINERS_FEM_TETRAHEDRON_BASIC_H_INCLUDED
