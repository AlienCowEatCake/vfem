#if !defined(CONTAINERS_FEM_TRIANGLE_BASIC_3D_H_INCLUDED)
#define CONTAINERS_FEM_TRIANGLE_BASIC_3D_H_INCLUDED

#include <cmath>
#include "triangle_basic.h"
#include "../generic/array_t.h"
#include "../generic/matrix_t.h"
#include "../geometry/point3_t.h"
#include "../geometry/vector3_t.h"
#include "../../cubatures/triangle_integration.h"

namespace core { namespace containers { namespace fem {

/**
 * @brief Класс треугольник в трехмерном пространстве
 */
template<typename edge, typename face, typename phys_area>
class triangle_basic_3d : public virtual triangle_basic<geometry::point3_t<double>, edge, face, phys_area>
{
public:
    triangle_basic_3d()
    {
        m_jacobian = 0.0;
    }

    /**
     * @brief Инициализация для работы в трехмерном пространстве
     * @param[in] tr_integration Параметры интегрирования
     */
    void init_3d(const cubatures::triangle_integration & tr_integration)
    {
        // Построение локальной системы координат
        geometry::vector3_t<double> g1(this->get_node(0), this->get_node(1));
        geometry::vector3_t<double> g2(this->get_node(0), this->get_node(2));
        geometry::vector3_t<double> e1 = g1 / g1.norm();
        geometry::vector3_t<double> e2 = g2 - ((g1 * g2) / (g1 * g1)) * g1;
        e2 = e2 / e2.norm();
        geometry::vector3_t<double> e3 = g1.cross(e2);
        e3 = e3 / e3.norm();

        for(std::size_t i = 0; i < 3; i++)
        {
            m_transition_matrix[0][i] = e1[i];
            m_transition_matrix[1][i] = e2[i];
            m_transition_matrix[2][i] = e3[i];
        }

        geometry::point3_t<double> local_nodes[3];
        local_nodes[0] = geometry::point3_t<double>(0.0, 0.0, 0.0);
        local_nodes[1] = geometry::point3_t<double>(g1.norm(), 0.0, 0.0);
        local_nodes[2] = (m_transition_matrix * g2).to_point();

        generic::matrix_t<double, 3, 3> D;
        // Формирование L-координат
        for(std::size_t i = 0; i < 2; i++)
            for(std::size_t j = 0; j < 3; j++)
                D[i][j] = local_nodes[j][i];

        for(std::size_t i = 0; i < 3; i++)
            D[2][i] = 1.0;

        double D_det;
        m_L = inverse(D, D_det);

        m_jacobian = fabs(D_det);

        // Точки Гаусса в локальной системе координат
        generic::array_t<geometry::point3_t<double> > gauss_points_local(tr_integration.get_gauss_num());
        for(std::size_t j = 0 ; j < tr_integration.get_gauss_num(); j++)
        {
            for(std::size_t i = 0; i < 2; i++)
            {
                gauss_points_local[j][i] = 0.0;
                for(size_t k = 0; k < 3; k++)
                    gauss_points_local[j][i] += D[i][k] * tr_integration.get_gauss_point_master(j, k);
            }
            gauss_points_local[j][2] = 0.0;
        }

        // Точки Гаусса в глобальной системе координат
        m_gauss_points.resize(tr_integration.get_gauss_num());
        for(std::size_t i = 0; i < tr_integration.get_gauss_num(); i++)
            m_gauss_points[i] = to_global(gauss_points_local[i]);
    }

    /**
     * @brief Получить значение L-координаты
     * @param[in] i Номер L-координаты
     * @param[in] p Точка, в которой нужно узнать значение
     * @return Значение L-координаты
     */
    double lambda(std::size_t i, const geometry::point3_t<double> & p) const
    {
        double result = 0.0;
        for(size_t j = 0; j < 2; j++)
            result += m_L[i][j] * p[j];
        result += m_L[i][2];
        return result;
    }

    /**
     * @brief Получить градиент L-координаты (в глобальных координатах)
     * @param[in] i Номер L-координаты
     * @return Градиент L-координаты в глобальных координатах
     */
    geometry::vector3_t<double> grad_lambda(std::size_t i) const
    {
        geometry::vector3_t<double> grad(m_L[i][0], m_L[i][1], 0.0);
        return to_global(grad);
    }

    /**
     * @brief Перевод точки в локальную систему координат треугольника
     * @param[in] p Точка в глобальной системе координат
     * @return Точка в локальной системе координат
     */
    geometry::point3_t<double> to_local(const geometry::point3_t<double> & p) const
    {
        geometry::point3_t<double> shift = p;
        // Cдвиг
        for(std::size_t i = 0; i < 3; i++)
            shift[i] -= this->get_node(0)[i];
        geometry::point3_t<double> turn;
        // Поворот
        for(std::size_t i = 0; i < 3; i++)
        {
            turn[i] = 0.0;
            for(std::size_t j = 0; j < 3; j++)
                turn[i] += m_transition_matrix[i][j] * shift[j];
        }
        return turn;
    }

    /**
     * @brief Перевод точки в глобальную систему координат
     * @param[in] p Точка в локальной системе координат
     * @return Точка в глобальной системе координат
     */
    geometry::point3_t<double> to_global(const geometry::point3_t<double> & p) const
    {
        geometry::point3_t<double> turn;
        // Поворот
        for(std::size_t i = 0; i < 3; i++)
        {
            turn[i] = 0.0;
            for(std::size_t j = 0; j < 3; j++)
                turn[i] += m_transition_matrix[j][i] * p[j];
        }
        geometry::point3_t<double> shift;
        // Сдвиг
        for(std::size_t i = 0; i < 3; i++)
            shift[i] = turn[i] + this->get_node(0)[i];
        return shift;
    }

    /**
     * @brief Перевод вектора в глобальную систему координат
     * @param[in] p Вектор в локальной системе координат
     * @return Вектор в глобальной системе координат
     */
    geometry::vector3_t<double> to_global(const geometry::vector3_t<double> & v) const
    {
        geometry::vector3_t<double> result;
        // Поворот
        for(std::size_t i = 0; i < 3; i++)
        {
            result[i] = 0;
            for(std::size_t j = 0; j < 3; j++)
                result[i] += m_transition_matrix[j][i] * v[j];
        }
        return result;
    }

    /**
     * @brief Получить якобиан для интегрирования
     * @return Якобиан
     */
    inline double get_jacobian() const
    {
        return m_jacobian;
    }

    /**
     * @brief Получить точку интегрирования (точку Гаусса)
     * @param[in] i Номер точки интегрирования
     * @return Точка интегрирования
     */
    inline const geometry::point3_t<double> & get_gauss_point(std::size_t i) const
    {
        return m_gauss_points[i];
    }

private:

    /**
     * @brief Матрица L-координат
     */
    generic::matrix_t<double> m_L;

    /**
     * @brief Матрица перехода между локальной и глобальной с.к.
     */
    generic::matrix_t<double, 3, 3> m_transition_matrix;

    /**
     * @brief Точки Гаусса
     */
    generic::array_t<geometry::point3_t<double> > m_gauss_points;

    /**
     * @brief Якобиан
     */
    double m_jacobian;
};

}}} // namespace core::containers::fem

#endif // CONTAINERS_FEM_TRIANGLE_BASIC_3D_H_INCLUDED
