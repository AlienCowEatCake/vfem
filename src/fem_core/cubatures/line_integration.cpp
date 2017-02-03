#include "line_integration.h"

#include <cmath>

namespace fem_core { namespace cubatures {

// *************************************************************************************************

namespace {

// Численное интегрирование на отрезках
// https://ru.wikipedia.org/wiki/Список_квадратурных_формул

namespace line_integration_3
{
    const std::size_t gauss_num = 2;
    const double gauss_weights[gauss_num] =
    {
        1.0,
        1.0
    };
    const double gauss_points_master[gauss_num] =
    {
        - 1.0 / std::sqrt(3.0),
        1.0 / std::sqrt(3.0)
    };
}

namespace line_integration_5
{
    const std::size_t gauss_num = 3;
    const double gauss_weights[gauss_num] =
    {
        5.0 / 9.0,
        8.0 / 9.0,
        5.0 / 9.0
    };
    const double gauss_points_master[gauss_num] =
    {
        - std::sqrt(3.0 / 5.0),
        0.0,
        std::sqrt(3.0 / 5.0)
    };
}

namespace line_integration_7
{
    const std::size_t gauss_num = 4;
    const double gauss_weights[gauss_num] =
    {
        (18.0 + std::sqrt(30.0)) / 36.0,
        (18.0 + std::sqrt(30.0)) / 36.0,
        (18.0 - std::sqrt(30.0)) / 36.0,
        (18.0 - std::sqrt(30.0)) / 36.0
    };
    const double gauss_points_master[gauss_num] =
    {
        - std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
        std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
        - std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
        std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0))
    };
}

namespace line_integration_9
{
    const std::size_t gauss_num = 5;
    const double gauss_weights[gauss_num] =
    {
        128.0 / 225.0,
        (322.0 + 13.0 * sqrt(70.0)) / 900.0,
        (322.0 + 13.0 * sqrt(70.0)) / 900.0,
        (322.0 - 13.0 * sqrt(70.0)) / 900.0,
        (322.0 - 13.0 * sqrt(70.0)) / 900.0
    };
    const double gauss_points_master[gauss_num] =
    {
        0.0,
        sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0,
        -sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0,
        sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0,
        -sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0
    };
}

} // namespace

// *************************************************************************************************

void line_integration::set(std::size_t num, const double weights[], const double points[], std::size_t num_steps)
{
    m_gauss_num = num * num_steps;
    m_gauss_weights.resize(m_gauss_num);
    m_gauss_points_master.resize(m_gauss_num);

    if(num_steps < 2)
    {
        for(std::size_t i = 0; i < m_gauss_num; i++)
        {
            m_gauss_weights[i] = weights[i] / 2.0;
            m_gauss_points_master[i] = points[i];
        }
    }
    else
    {
        const double h = 2.0 / static_cast<double>(num_steps);
        for(std::size_t k = 0, i = 0; k < num_steps; k++)
        {
            const double begin = h * static_cast<double>(k);
            for(std::size_t j = 0; j < num; j++)
            {
                m_gauss_weights[i] = weights[j] * h / 4.0;
                m_gauss_points_master[i] = (points[j] + 1.0) / 2.0 * h + begin - 1.0;
                i++;
            }
        }
    }
}

line_integration::line_integration(std::size_t order, std::size_t num_steps)
{
    init(order, num_steps);
}

void line_integration::init(std::size_t order, std::size_t num_steps)
{
    if(order < 3)
        order = 3;
    if(order > 9)
        order = 9;

    switch(order)
    {
    case 3:
    {
        namespace curr = line_integration_3;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master, num_steps);
        break;
    }
    case 4:
    case 5:
    {
        namespace curr = line_integration_5;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master, num_steps);
        break;
    }
    case 6:
    case 7:
    {
        namespace curr = line_integration_7;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master, num_steps);
        break;
    }
    case 8:
    case 9:
    {
        namespace curr = line_integration_9;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master, num_steps);
        break;
    }
    default:
    {
        namespace curr = line_integration_9;
        set(curr::gauss_num, curr::gauss_weights, curr::gauss_points_master, num_steps);
        break;
    }
    }
}

// *************************************************************************************************

}} // namespace fem_core::cubatures
