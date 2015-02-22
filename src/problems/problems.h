#if !defined PROBLEMS_H_INCLUDED
#define PROBLEMS_H_INCLUDED

#include "../vfem/vfem.h"

//#define ANALYTICAL_CUBE
//#define FIELD_CUBE
//#define AREA_POINT_SOURCE
//#define AREA_SOURCE
//#define CUBE_PML
//#define AREA_PML
//#define SOURCE_PML
#define AREA_PML_SOURCE

extern string mesh_filename;
extern string phys_filename;
extern string tecplot_filename;

void postprocessing(VFEM & v, char * timebuf);

#endif // PROBLEMS_H_INCLUDED
