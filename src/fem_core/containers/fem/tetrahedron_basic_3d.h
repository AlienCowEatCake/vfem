#if !defined(CONTAINERS_FEM_TETRAHEDRON_BASIC_3D_H_INCLUDED)
#define CONTAINERS_FEM_TETRAHEDRON_BASIC_3D_H_INCLUDED

#include <cmath>
#include "tetrahedron_basic.h"
#include "../generic/array_t.h"
#include "../generic/matrix_t.h"
#include "../geometry/point3_t.h"
#include "../geometry/vector3_t.h"
#include "../../cubatures/tetrahedron_integration.h"

namespace fem_core { namespace containers { namespace fem {

/**
 * @brief Класс тетраэдр в трехмерном пространстве
 */
template<typename edge, typename face, typename phys_area>
class tetrahedron_basic_3d : public tetrahedron_basic<geometry::point3_t<double>, edge, face, phys_area>
{
public:
    tetrahedron_basic_3d()
        : m_jacobian(0.0)
    {}

    /**
     * @brief Инициализация для работы в трехмерном пространстве
     * @param[in] tet_integration Параметры интегрирования
     */
    void init_3d(const cubatures::tetrahedron_integration & tet_integration)
    {
        generic::matrix_t<double, 4, 4> D;
        for(std::size_t i = 0; i < 3; i++)
            for(std::size_t j = 0; j < 4; j++)
                D[i][j] = this->get_node(j)[i];
        for(std::size_t j = 0; j < 4; j++)
            D[3][j] = 1.0;
        double D_det;
        m_L = inverse(D, D_det);

        m_jacobian = std::fabs(D_det);

        // Перевод точек Гаусса с мастер-элемента на текущий тетраэдр
        m_gauss_points.resize(tet_integration.get_gauss_num());
        for(std::size_t i = 0; i < 3; i++)
        {
            for(std::size_t j = 0 ; j < tet_integration.get_gauss_num(); j++)
            {
                m_gauss_points[j][i] = 0;
                for(std::size_t k = 0; k < 4; k++)
                    m_gauss_points[j][i] += D[i][k] * tet_integration.get_gauss_point_master(j, k);
            }
        }

        // Расчет барицентра
        for(std::size_t i = 0; i < 3; i++)
        {
            m_barycenter[i] = 0.0;
            for(std::size_t j = 0; j < 4; j++)
                m_barycenter[i] += this->get_node(j)[i];
            m_barycenter[i] /= 4.0;
        }

        // Представление прямых тетраэдра в параметрическом виде a*t + b
        // Само ребо получается при 0<=t<=1
        m_edges_a.resize(6, 3);
        m_edges_b.resize(6, 3);
        for(std::size_t i = 0; i < 3; i++)
        {
            m_edges_a[0][i] = this->get_node(1)[i] - this->get_node(0)[i];
            m_edges_b[0][i] = this->get_node(0)[i];
            m_edges_a[1][i] = this->get_node(2)[i] - this->get_node(0)[i];
            m_edges_b[1][i] = this->get_node(0)[i];
            m_edges_a[2][i] = this->get_node(3)[i] - this->get_node(0)[i];
            m_edges_b[2][i] = this->get_node(0)[i];
            m_edges_a[3][i] = this->get_node(2)[i] - this->get_node(1)[i];
            m_edges_b[3][i] = this->get_node(1)[i];
            m_edges_a[4][i] = this->get_node(3)[i] - this->get_node(1)[i];
            m_edges_b[4][i] = this->get_node(1)[i];
            m_edges_a[5][i] = this->get_node(3)[i] - this->get_node(2)[i];
            m_edges_b[5][i] = this->get_node(2)[i];
        }
    }

    /**
     * @brief Определить, внутри тетраэдра точка или нет
     * @param[in] p Точка
     * @return true, если внутри ; false, если иначе
     */
    bool inside(const geometry::point3_t<double> & p) const
    {
        for(std::size_t i = 0; i < 4; i++)
        {
            const double a = lambda(i, p);
            if(a -1e-12 > 1.0 || a +1e-12 < 0.0)
                return false;
        }
        return true;
    }

    /**
     * @brief Определить, внутри тетраэдра точка или нет
     * @param[in] x Координата x точки
     * @param[in] y Координата y точки
     * @param[in] z Координата z точки
     * @return true, если внутри ; false, если иначе
     */
    bool inside(double x, double y, double z) const
    {
        return inside(geometry::point3_t<double>(x, y, z));
    }

    /**
     * @brief Определить, лежит ли тетраэдр внутри куба или нет
     * @param[in] x0 Координата x начала куба
     * @param[in] x1 Координата x конца куба
     * @param[in] y0 Координата y начала куба
     * @param[in] y1 Координата y конца куба
     * @param[in] z0 Координата z начала куба
     * @param[in] z1 Координата z конца куба
     * @return true, если внутри ; false, если иначе
     */
    bool in_cube(double x0, double x1, double y0, double y1, double z0, double z1) const
    {
        // Тетраэдр внутри куба
        if(get_barycenter().inside(x0, x1, y0, y1, z0, z1))
            return true;
        for(std::size_t i = 0; i < 4; i++)
            if(this->get_node(i).inside(x0, x1, y0, y1, z0, z1))
                return true;

        // Куб внутри тетраэдра
        if(inside(x0, y0, z0) || inside(x1, y0, z0) || inside(x0, y1, z0) || inside(x1, y1, z0) ||
                inside(x0, y0, z1) || inside(x1, y0, z1) || inside(x0, y1, z1) || inside(x1, y1, z1))
            return true;

        // Пересечение куба и тетраэдра, не обработанное выше
        double t[6][6];             // Параметры прямых для всех прямых (6) и всех плоскостей (6)
        double cords_t[6][6][3];    // Прямая, плоскость, координата
        for(std::size_t i = 0; i < 6; i++)
        {
            t[i][0] = (x0 - m_edges_b[i][0]) / m_edges_a[i][0]; // x = x0
            t[i][1] = (x1 - m_edges_b[i][0]) / m_edges_a[i][0]; // x = x1
            t[i][2] = (y0 - m_edges_b[i][1]) / m_edges_a[i][1]; // y = y0
            t[i][3] = (y1 - m_edges_b[i][1]) / m_edges_a[i][1]; // y = y1
            t[i][4] = (z0 - m_edges_b[i][2]) / m_edges_a[i][2]; // z = z0
            t[i][5] = (z1 - m_edges_b[i][2]) / m_edges_a[i][2]; // z = z1
        }

        for(std::size_t i = 0; i < 6; i++)
            for(std::size_t j = 0; j < 6; j++)
                for(std::size_t k = 0; k < 3; k++)
                    cords_t[i][j][k] = m_edges_a[i][k] * t[i][j] + m_edges_b[i][k];

        for(std::size_t i = 0; i < 6; i++)  // Берем прямую и проверяем, что пересечение с плоскостью попадает в раcсматриваемый отрезок прямой и в рассатриваемую часть плоскости
        {
            if(
                (t[i][0] >= 0.0 && t[i][0] <= 1.0 && cords_t[i][0][1] >= y0 && cords_t[i][0][1] >= y1 && cords_t[i][0][2] >= z0 && cords_t[i][0][2] <= z1) || // x = x0
                (t[i][1] >= 0.0 && t[i][1] <= 1.0 && cords_t[i][1][1] >= y0 && cords_t[i][1][1] >= y1 && cords_t[i][1][2] >= z0 && cords_t[i][1][2] <= z1) || // x = x1
                (t[i][2] >= 0.0 && t[i][2] <= 1.0 && cords_t[i][2][0] >= x0 && cords_t[i][2][0] <= x1 && cords_t[i][2][2] >= z0 && cords_t[i][2][2] <= z1) || // y = y0
                (t[i][3] >= 0.0 && t[i][3] <= 1.0 && cords_t[i][3][0] >= x0 && cords_t[i][3][0] <= x1 && cords_t[i][3][2] >= z0 && cords_t[i][3][2] <= z1) || // y = y1
                (t[i][4] >= 0.0 && t[i][4] <= 1.0 && cords_t[i][4][0] >= x0 && cords_t[i][4][0] <= x1 && cords_t[i][4][1] >= y0 && cords_t[i][4][1] <= y1) || // z = z0
                (t[i][5] >= 0.0 && t[i][5] <= 1.0 && cords_t[i][5][0] >= x0 && cords_t[i][5][0] <= x1 && cords_t[i][5][1] >= y0 && cords_t[i][5][1] <= y1)    // z = z1
            )
                return true;
        }
        return false;
    }

    /**
     * @brief Получить барицентр тетраэдра
     * @return Барицентр тетраэдра
     */
    inline const geometry::point3_t<double> & get_barycenter() const
    {
        return m_barycenter;
    }

    /**
     * @brief Получить градиент L-координаты
     * @param[in] i Номер L-координаты
     * @return Градиент L-координаты
     */
    geometry::vector3_t<double> grad_lambda(std::size_t i) const
    {
        assert(i < 4); // Если i == 4, то явно где-то косяк
        return geometry::vector3_t<double>(m_L[i][0], m_L[i][1], m_L[i][2]);
    }

    /**
     * @brief Получить значение L-координаты
     * @param[in] i Номер L-координаты
     * @param[in] p Точка, в которой нужно узнать значение
     * @return Значение L-координаты
     */
    double lambda(std::size_t i, const geometry::point3_t<double> & p) const
    {
        assert(i < 4);
        return m_L[i][3] + m_L[i][0] * p.x + m_L[i][1] * p.y + m_L[i][2] * p.z;
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
     * @brief Барицентр тетраэдра
     */
    geometry::point3_t<double> m_barycenter;

    /**
     * @brief Матрица L-координат
     */
    generic::matrix_t<double> m_L;

    /**
     * @brief Точки Гаусса
     */
    generic::array_t<geometry::point3_t<double> > m_gauss_points;

    /**
     * @brief Якобиан
     */
    double m_jacobian;

    /**
     * @brief Параметры прямых для дерева
     */
    generic::matrix_t<double> m_edges_a, m_edges_b;
};

}}} // namespace fem_core::containers::fem

#endif // CONTAINERS_FEM_TETRAHEDRON_BASIC_3D_H_INCLUDED
