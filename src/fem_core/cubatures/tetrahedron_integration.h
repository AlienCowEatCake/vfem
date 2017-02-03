#if !defined(CUBATURES_TETRAHEDRON_INTEGRATION_H_INCLUDED)
#define CUBATURES_TETRAHEDRON_INTEGRATION_H_INCLUDED

#include <cstddef>
#include "../containers/generic/array_t.h"
#include "../containers/generic/matrix_t.h"

namespace fem_core { namespace cubatures {

/**
 * @brief Класс для хранения точек интегрирования на тетраэдрах
 * @note Мастер-элемент - тетраэдр с вершинами (0,0,0), (1,0,0), (0,1,0) и (0,0,1)
 */
class tetrahedron_integration
{
public:
    tetrahedron_integration(std::size_t order = 8);

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
     * @note Веса заданы с учетом объема мастер-элемента, в вычислении якобиана его учитывать не требуется
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
        assert(coord_num < 4);
        return m_gauss_points_master[point_num][coord_num];
    }

    /**
     * @brief Инициализация структуры для заданного порядка интегрирования
     * @param[in] order Порядок интегрирования (от 2 до 14 включительно)
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
     * @brief Точки интегрирования на мастер-тетраэдре, в каждой строке по 4 L-координаты
     */
    containers::generic::matrix_t<double> m_gauss_points_master;

    /**
     * @brief Конвертер для внутренних нужд
     */
    void set(std::size_t num, const double weights[], const double points[][4]);
};

}} // namespace fem_core::cubatures

#endif // CUBATURES_TETRAHEDRON_INTEGRATION_H_INCLUDED
