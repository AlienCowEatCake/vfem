#include "problems.h"

// Простое сравнение (на геометрически идентичных подобластях сеток и одинаковых базисах)
void compare_simple(VFEM & master, VFEM & slave, vector<diff_area> & areas)
{
    double diff = 0.0, norm = 0.0;
    vector3 diff_v3(0.0, 0.0, 0.0), norm_v3(0.0, 0.0, 0.0);
    cvector3 diff_cv3(0.0, 0.0, 0.0), norm_cv3(0.0, 0.0, 0.0);
    for(size_t i = 0; i < master.fes.size(); i++)
    {
        // Флаг обработки
        bool flag = false;
        // Поймем, нужно ли обрабатывать текущий КЭ
        finite_element * fe_m = &(master.fes[i]);
        for(size_t j = 0; j < areas.size(); j++)
        {
            if(fe_m->in_cube(areas[j].p1.x, areas[j].p2.x, areas[j].p1.y, areas[j].p2.y, areas[j].p1.z, areas[j].p2.z))
            {
                if(areas[j].included)
                    flag = true;
                else
                    break;
            }
        }
        // Если нужно - обработаем
        if(flag)
        {
            finite_element * fe_s = slave.get_fe(fe_m->barycenter);
            assert(fe_s);

            // Находим локальные веса
            array_t<complex<double> > q_m(master.config.basis.tet_bf_num);
            array_t<complex<double> > q_s(master.config.basis.tet_bf_num);
            for(size_t k = 0; k < master.config.basis.tet_bf_num; k++)
            {
                q_m[i] = master.slae.x[master.get_tet_dof(fe_m, k)];
                q_s[i] = slave.slae.x[slave.get_tet_dof(fe_s, k)];
            }

            // Считаем разность в норме L2
            trio<double, vector3, cvector3> diff_tmp = fe_m->diff_normL2(q_m, q_s);
            trio<double, vector3, cvector3> norm_tmp = fe_m->normL2(q_s);
            diff += diff_tmp.first;
            diff_v3 += diff_tmp.second;
            diff_cv3 += diff_tmp.third;
            norm += norm_tmp.first;
            norm_v3 += norm_tmp.second;
            norm_cv3 += norm_tmp.third;
        }
    }
    cout << "Diff (L2):        " << sqrt(diff / norm) << endl;
    cout << "Diff (L2) [Ex]:   " << sqrt(diff_v3.x / norm_v3.x) << endl;
    cout << "Diff (L2) [Ey]:   " << sqrt(diff_v3.y / norm_v3.y) << endl;
    cout << "Diff (L2) [Ez]:   " << sqrt(diff_v3.z / norm_v3.z) << endl;
    cout << "Diff (L2) [ExR]:  " << sqrt(diff_cv3.x.real() / norm_cv3.x.real()) << endl;
    cout << "Diff (L2) [EyR]:  " << sqrt(diff_cv3.y.real() / norm_cv3.y.real()) << endl;
    cout << "Diff (L2) [EzR]:  " << sqrt(diff_cv3.z.real() / norm_cv3.z.real()) << endl;
    cout << "Diff (L2) [ExI]:  " << sqrt(diff_cv3.x.imag() / norm_cv3.x.imag()) << endl;
    cout << "Diff (L2) [EyI]:  " << sqrt(diff_cv3.y.imag() / norm_cv3.y.imag()) << endl;
    cout << "Diff (L2) [EzI]:  " << sqrt(diff_cv3.z.imag() / norm_cv3.z.imag()) << endl;
}

// Функция для получения решения (прототипом соответствующая func_true) для сравнения решений
cvector3 func_diff(const point & p, const phys_area & phys, void * data)
{
    MAYBE_UNUSED(phys);
    VFEM * v = (VFEM *)(data);
    return v->solution(p);
}

// Сложное сравнение (на произвольных подобластях сеток и базисах)
void compare_complex(VFEM & master, VFEM & slave, vector<diff_area> & areas)
{
    double diff = 0.0, norm = 0.0;
    vector3 diff_v3(0.0, 0.0, 0.0), norm_v3(0.0, 0.0, 0.0);
    cvector3 diff_cv3(0.0, 0.0, 0.0), norm_cv3(0.0, 0.0, 0.0);
    for(size_t i = 0; i < master.fes.size(); i++)
    {
        // Флаг обработки
        bool flag = false;
        // Поймем, нужно ли обрабатывать текущий КЭ
        finite_element * fe_m = &(master.fes[i]);
        for(size_t j = 0; j < areas.size(); j++)
        {
            if(fe_m->in_cube(areas[j].p1.x, areas[j].p2.x, areas[j].p1.y, areas[j].p2.y, areas[j].p1.z, areas[j].p2.z))
            {
                if(areas[j].included)
                    flag = true;
                else
                    break;
            }
        }
        // Если нужно - обработаем
        if(flag)
        {
            // Находим локальные веса
            array_t<complex<double> > q_m(master.config.basis.tet_bf_num);
            for(size_t k = 0; k < master.config.basis.tet_bf_num; k++)
                q_m[i] = master.slae.x[master.get_tet_dof(fe_m, k)];

            // Считаем разность в норме L2
            trio<double, vector3, cvector3> diff_tmp = fe_m->diff_normL2(q_m, func_diff, & slave);
            trio<double, vector3, cvector3> norm_tmp = fe_m->normL2(func_diff, & slave);
            diff += diff_tmp.first;
            diff_v3 += diff_tmp.second;
            diff_cv3 += diff_tmp.third;
            norm += norm_tmp.first;
            norm_v3 += norm_tmp.second;
            norm_cv3 += norm_tmp.third;
        }
    }
    cout << "Diff (L2):        " << sqrt(diff / norm) << endl;
    cout << "Diff (L2) [Ex]:   " << sqrt(diff_v3.x / norm_v3.x) << endl;
    cout << "Diff (L2) [Ey]:   " << sqrt(diff_v3.y / norm_v3.y) << endl;
    cout << "Diff (L2) [Ez]:   " << sqrt(diff_v3.z / norm_v3.z) << endl;
    cout << "Diff (L2) [ExR]:  " << sqrt(diff_cv3.x.real() / norm_cv3.x.real()) << endl;
    cout << "Diff (L2) [EyR]:  " << sqrt(diff_cv3.y.real() / norm_cv3.y.real()) << endl;
    cout << "Diff (L2) [EzR]:  " << sqrt(diff_cv3.z.real() / norm_cv3.z.real()) << endl;
    cout << "Diff (L2) [ExI]:  " << sqrt(diff_cv3.x.imag() / norm_cv3.x.imag()) << endl;
    cout << "Diff (L2) [EyI]:  " << sqrt(diff_cv3.y.imag() / norm_cv3.y.imag()) << endl;
    cout << "Diff (L2) [EzI]:  " << sqrt(diff_cv3.z.imag() / norm_cv3.z.imag()) << endl;
}
