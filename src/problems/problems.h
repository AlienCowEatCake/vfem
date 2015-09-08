#if !defined PROBLEMS_H_INCLUDED
#define PROBLEMS_H_INCLUDED

#include "../vfem/vfem.h"

//#define ANALYTICAL_CUBE
//#define LOOP_PML
//#define AREA_2LAYERS_LOOP_PML
//#define AREA_2LAYERS_LOOP_MANY_PML
//#define AREA_2LAYERS_LOOP_UNIVERSAL_PML
//#define AREA_3LAYERS_INC_LOOP_PML
//#define AREA_4LAYERS_LOOP_PML
#define PROBLEM_STANDARD
//#define PROBLEM_STANDARD_PML

void postprocessing(VFEM & v, char * timebuf);

#endif // PROBLEMS_H_INCLUDED
