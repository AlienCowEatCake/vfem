#if !defined(CUBATURES_LINE_INTEGRATION_H_INCLUDED)
#define CUBATURES_LINE_INTEGRATION_H_INCLUDED

#include <cstddef>
#include "../containers/generic/array_t.h"

namespace fem_core { namespace cubatures {

/**
 * @brief Класс для хранения точек интегрирования на отрезках
 * @note Мастер-элемент - отрезок [-1,1]
 */
class line_integration
{
public:
    /**
     * @brief Конструктор
     * @param[in] order Порядок интегрирования (от 3 до 9 включительно)
     * @param[in] num_steps На сколько отрезков требуется разбить исходный
     */
    line_integration(std::size_t order = 9, std::size_t num_steps = 1);

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
     * @note Веса заданы с учетом длины мастер-элемента, в вычислении якобиана его учитывать не требуется
     */
    inline double get_gauss_weight(std::size_t point_num) const
    {
        assert(point_num < m_gauss_num);
        return m_gauss_weights[point_num];
    }

    /**
     * @brief Получить координату точки интегрирования
     * @param[in] point_num Номер точки интегрирования
     * @return Координата искомой точки интегрирования
     */
    inline double get_gauss_point_master(std::size_t point_num) const
    {
        assert(point_num < m_gauss_num);
        return m_gauss_points_master[point_num];
    }

    /**
     * @brief Инициализация структуры для заданного порядка интегрирования
     * @param[in] order Порядок интегрирования (от 3 до 9 включительно)
     * @param[in] num_steps На сколько отрезков требуется разбить исходный
     */
    void init(std::size_t order, std::size_t num_steps = 1);

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
     * @brief Точки интегрирования на мастер-отрезке
     */
    containers::generic::array_t<double> m_gauss_points_master;

    /**
     * @brief Конвертер для внутренних нужд
     */
    void set(std::size_t num, const double weights[], const double points[], std::size_t num_steps);
};

}} // namespace fem_core::cubatures

#endif // CUBATURES_LINE_INTEGRATION_H_INCLUDED
