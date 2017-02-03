#if !defined(CUBATURES_TRIANGLE_INTEGRATION_H_INCLUDED)
#define CUBATURES_TRIANGLE_INTEGRATION_H_INCLUDED

#include <cstddef>
#include <cassert>
#include "../containers/generic/array_t.h"
#include "../containers/generic/matrix_t.h"

namespace fem_core { namespace cubatures {

/**
 * @brief Класс для хранения точек интегрирования на треугольниках
 * @note Мастер-элемент - треугольник с вершинами (0,0), (1,0) и (0,1)
 */
class triangle_integration
{
public:
    triangle_integration(std::size_t order = 8);

    /**
     * @brief Получить количество точек интегрирования
     * @return Количество точек интегрирования
     */
    inline std::size_t get_gauss_num() const
    {
        return m_gauss_num;
    }

    /**
     * @brief Получить вес одной из точек интегрирования
     * @param[in] point_num
     * @return
     * @note Веса заданы с учетом площади мастер-элемента, в вычислении якобиана ее учитывать не требуется
     */
    inline double get_gauss_weight(std::size_t point_num) const
    {
        assert(point_num < m_gauss_num);
        return m_gauss_weights[point_num];
    }

    /**
     * @brief Получить одну из координат одной из точек интегрирования
     * @param[in] point_num Номер точки интегрирования
     * @param[in] coord_num Номер координаты интегрирования
     * @return Координата искомой точки интегрирования
     */
    inline double get_gauss_point_master(std::size_t point_num, std::size_t coord_num) const
    {
        assert(point_num < m_gauss_num);
        assert(coord_num < 3);
        return m_gauss_points_master[point_num][coord_num];
    }

    /**
     * @brief Инициализация структуры для заданного порядка интегрирования
     * @param[in] order Порядок интегрирования (от 2 до 29 включительно)
     */
    void init(std::size_t order);

private:

    /**
     * @brief Количество точек интегрирования
     */
    std::size_t m_gauss_num;

    /**
     * @brief Веса интегрирования
     */
    containers::generic::array_t<double> m_gauss_weights;

    /**
     * @brief Точки интегрирования на мастер-треугольнике, в каждой строке по 3 L-координаты
     */
    containers::generic::matrix_t<double> m_gauss_points_master;

    /**
     * @brief Конвертер для внутренних нужд
     */
    void set(std::size_t num, const double weights[], const double points[][3]);
};

}} // namespace fem_core::cubatures

#endif // CUBATURES_TRIANGLE_INTEGRATION_H_INCLUDED
