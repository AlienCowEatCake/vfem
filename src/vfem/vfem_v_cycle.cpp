#include "vfem.h"

// Проектирование на пространство ядра
void VFEM::to_kernel_space(const complex<double> * in, complex<double> * out) const
{
    size_t nodes_num = nodes.size();
    size_t edges_num = edges.size();
#if BASIS_ORDER >= 2
    size_t faces_num = faces.size();
#endif

    for(size_t i = 0; i < nodes_num; i++)
        out[i] = 0.0;

    for(set<edge>::iterator it = edges.begin(); it != edges.end(); ++it)
    {
        out[it->nodes[0]->num] -= in[it->num];
        out[it->nodes[1]->num] += in[it->num];
    }

#if BASIS_ORDER >= 2 || BASIS_TYPE == 2
    for(set<edge>::iterator it = edges.begin(); it != edges.end(); ++it)
        out[nodes_num + it->num] = in[edges_num + it->num];
#endif

#if BASIS_ORDER > 2 || (BASIS_TYPE == 2 && BASIS_ORDER == 2)
    for(set<face>::iterator it = faces.begin(); it != faces.end(); ++it)
        out[nodes_num + edges_num + it->num] = in[2 * edges_num + 2 * faces_num + it->num];
    for(set<edge>::iterator it = edges.begin(); it != edges.end(); ++it)
        out[nodes_num + edges_num + faces_num + it->num] = in[2 * edges_num + 3 * faces_num + it->num];
#endif
}

// Интерполяция на полное пространство
void VFEM::to_full_space(const complex<double> * in, complex<double> * out) const
{
    size_t nodes_num = nodes.size();
    size_t edges_num = edges.size();
#if BASIS_ORDER >= 2
    size_t faces_num = faces.size();
#endif

    for(size_t i = 0; i < edges_num; i++)
        out[i] = 0.0;

    for(set<edge>::iterator it = edges.begin(); it != edges.end(); ++it)
    {
        out[it->num] -= in[it->nodes[0]->num];
        out[it->num] += in[it->nodes[1]->num];
    }

#if BASIS_ORDER >= 2 || BASIS_TYPE == 2
    for(set<edge>::iterator it = edges.begin(); it != edges.end(); ++it)
        out[edges_num + it->num] = in[nodes_num + it->num];
#endif

#if BASIS_ORDER > 2 || (BASIS_TYPE == 2 && BASIS_ORDER == 2)
    for(set<face>::iterator it = faces.begin(); it != faces.end(); ++it)
        out[2 * edges_num + 2 * faces_num + it->num] = in[nodes_num + edges_num + it->num];
    for(set<edge>::iterator it = edges.begin(); it != edges.end(); ++it)
        out[2 * edges_num + 3 * faces_num + it->num] = in[nodes_num + edges_num + faces_num + it->num];
#endif
}

// Запуск решения СЛАУ
void VFEM::solve()
{
    extern double SLAE_MAIN_EPSILON;
    slae.solve(SLAE_MAIN_EPSILON);
}
