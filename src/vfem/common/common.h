#if !defined(COMMON_H_INCLUDED)
#define COMMON_H_INCLUDED

#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif

#if defined(_MSC_VER)
#if !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#else
#include <cmath>
#endif

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cfloat>
#include <algorithm>
#include <cstring>
#include <vector>
#include <set>
#include <string>
#include <complex>
#include <map>
#include <utility>
#include <cassert>
#include <cctype>
#include <list>
#include <ctime>

#include "../../fem_core/containers/fem/triangle_basic_3d.h"
#include "../../fem_core/containers/fem/tetrahedron_basic_3d.h"
#include "../../fem_core/containers/generic/trio.h"
#include "../../fem_core/containers/generic/array_t.h"
#include "../../fem_core/containers/generic/matrix_t.h"
#include "../../fem_core/containers/geometry/edge_basic.h"
#include "../../fem_core/containers/geometry/face_triangle_basic.h"
#include "../../fem_core/containers/geometry/point3_t.h"
#include "../../fem_core/containers/geometry/vector3_t.h"
#include "../../fem_core/containers/tree/quadtree.h"
#include "../../fem_core/containers/tree/octree.h"
#include "../../fem_core/cubatures/triangle_integration.h"
#include "../../fem_core/cubatures/tetrahedron_integration.h"
#include "../../fem_core/cubatures/line_integration.h"
#include "../../fem_core/utils/cxxversion.h"
#include "../../fem_core/utils/fpu.h"
#include "../../fem_core/utils/nosighup.h"
#include "../../fem_core/utils/strings.h"
#include "../../fem_core/utils/progress.h"
#include "../../fem_core/utils/timers.h"
#include "../../fem_core/utils/inifile.h"
#include "../../fem_core/utils/evaluable_inifile.h"
#include "../../fem_core/evaluator/evaluator.h"
#include "../../fem_core/evaluator/evaluator_xyz.h"
#include "../../fem_core/solvers/CSLR/solvers_factory.h"
#include "../../fem_core/wrappers/omp_wrapper.h"

using namespace std;
using namespace fem_core::containers::fem;
using namespace fem_core::containers::generic;
using namespace fem_core::containers::geometry;
using namespace fem_core::containers::tree;
using namespace fem_core::cubatures;
using namespace fem_core::utils;
using namespace fem_core::utils::nosighup;
using namespace fem_core::utils::strings;
using namespace fem_core::utils::progress;
using namespace fem_core::utils::timers;
using namespace fem_core::utils::fpu;
using namespace fem_core::solvers::CSLR::factory;
using namespace fem_core::solvers::CSLR::preconditioners;
using namespace fem_core::solvers::CSLR::symmetric;
using namespace fem_core::wrappers::omp;

typedef point3_t<double> point;
typedef point3_t< complex<double> > cpoint;
typedef vector3_t<double> vector3;
typedef vector3_t< complex<double> > cvector3;

typedef face_triangle_basic<point> face;

#endif // COMMON_H_INCLUDED
