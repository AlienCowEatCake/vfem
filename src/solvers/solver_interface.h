#if !defined SOLVER_INTERFACE_H
#define SOLVER_INTERFACE_H

#include <cstdlib>
#include <iostream>
#include <typeinfo>

// Абстрактный решатель
template<typename val_type, typename ind_type = size_t>
class solver_interface
{
public:

    // Инициализация симметричная
    virtual void init(const ind_type * gi, const ind_type * gj, const val_type * di,
                      const val_type * gg, ind_type n) = 0;

    // Запуск решения
    virtual void solve(val_type * solution, const val_type * rp,
                       double eps, ind_type max_iter) = 0;

    // Деструктор
    virtual ~solver_interface() {}

};

#endif // SOLVER_INTERFACE_H

