#include "vfem.h"

template<typename T>
T avg(const T & v1, const T & v2)
{
    return 0.5 * (v1 + v2);
}

template<typename T>
T jmp_T(const T & v1, const T & n1, const T & v2, const T & n2)
{
    return n1.cross(v1) + n2.cross(v2);
}

cvector3 flux_E(const cvector3 & E1, const cvector3 & n1,
                const cvector3 & E2, const cvector3 & n2,
                const dg_config * dg_cfg)
{
    return avg(E1, E2)
            + dg_cfg->b * jmp_T(E1, n1, E2, n2);
}

cvector3 flux_rotE(const cvector3 & E1, const cvector3 & mu1rotE1, const cvector3 & n1,
                   const cvector3 & E2, const cvector3 & mu1rotE2, const cvector3 & n2,
                   const dg_config * dg_cfg)
{
    return avg(mu1rotE1, mu1rotE2)
            - dg_cfg->a * jmp_T(E1, n1, E2, n2)
            + dg_cfg->b * jmp_T(mu1rotE1, n1, mu1rotE2, n2);
}

// Применение численных потоков в DG
void VFEM::applying_dg_fluxes()
{
    for(size_t i = 0; i < config.dg.size(); i++)
    {
        dg_config * dg_cfg = &(config.dg[i]);
        size_t master_surface = dg_cfg->master_surface;
        size_t slave_surface = dg_cfg->slave_surface;

        vector<triangle_base>::iterator m_it = trs_dg[master_surface].begin();
        for(; m_it != trs_dg[master_surface].end(); ++m_it)
        {
            point barycenter_master;
            for(size_t i = 0; i < 3; i++)
            {
                barycenter_master[i] = 0.0;
                for(size_t j = 0; j < 3; j++)
                    barycenter_master[i] += m_it->get_node(j)[i];
                barycenter_master[i] /= 3.0;
            }

            vector<triangle_base>::iterator s_it = trs_dg[slave_surface].begin();
            for(; s_it != trs_dg[slave_surface].end(); ++s_it)
            {
                point barycenter_slave;
                for(size_t i = 0; i < 3; i++)
                {
                    barycenter_slave[i] = 0.0;
                    for(size_t j = 0; j < 3; j++)
                        barycenter_slave[i] += m_it->get_node(j)[i];
                    barycenter_slave[i] /= 3.0;
                }

                if(vector3(barycenter_master, barycenter_slave).norm2() < 1e-6)
                {




                    cout << "Here!" << endl;
                }
            }
        }
    }
}
