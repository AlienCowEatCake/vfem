#if !defined PROBLEMS_H_INCLUDED
#define PROBLEMS_H_INCLUDED

#include "../vfem/vfem.h"

class diff_area
{
public:
    bool included;
    point p1, p2;
};

void postprocessing(VFEM & v, char * timebuf);

#endif // PROBLEMS_H_INCLUDED
