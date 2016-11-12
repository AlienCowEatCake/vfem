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

#include "../../core/containers/fem/triangle_basic_3d.h"
#include "../../core/containers/fem/tetrahedron_basic_3d.h"
#include "../../core/containers/generic/trio.h"
#include "../../core/containers/generic/array_t.h"
#include "../../core/containers/generic/matrix_t.h"
#include "../../core/containers/geometry/edge_basic.h"
#include "../../core/containers/geometry/face_triangle_basic.h"
#include "../../core/containers/geometry/point3_t.h"
#include "../../core/containers/geometry/vector3_t.h"
#include "../../core/containers/tree/quadtree.h"
#include "../../core/containers/tree/octree.h"
#include "../../core/cubatures/triangle_integration.h"
#include "../../core/cubatures/tetrahedron_integration.h"
#include "../../core/cubatures/line_integration.h"
#include "../../core/utils/cxxversion.h"
#include "../../core/utils/fpu.h"
#include "../../core/utils/nosighup.h"
#include "../../core/utils/strings.h"
#include "../../core/utils/progress.h"
#include "../../core/utils/timers.h"
#include "../../core/utils/inifile.h"
#include "../../core/utils/evaluable_inifile.h"
#include "../../core/evaluator/evaluator.h"
#include "../../core/evaluator/evaluator_xyz.h"
#include "../../core/solvers/CSLR/solvers_factory.h"
#include "../../core/wrappers/omp_wrapper.h"

using namespace std;
using namespace core::containers::fem;
using namespace core::containers::generic;
using namespace core::containers::geometry;
using namespace core::containers::tree;
using namespace core::cubatures;
using namespace core::utils;
using namespace core::utils::nosighup;
using namespace core::utils::strings;
using namespace core::utils::progress;
using namespace core::utils::timers;
using namespace core::utils::fpu;
using namespace core::solvers::CSLR::factory;
using namespace core::solvers::CSLR::preconditioners;
using namespace core::solvers::CSLR::symmetric;
using namespace core::wrappers::omp;

typedef point3_t<double> point;
typedef point3_t< complex<double> > cpoint;
typedef vector3_t<double> vector3;
typedef vector3_t< complex<double> > cvector3;

typedef face_triangle_basic<point> face;

#endif // COMMON_H_INCLUDED
