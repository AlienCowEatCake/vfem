#if defined _MSC_VER && !defined _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "COCG_Di_Smooth_MKL.h"

#define COCG_LLT_Smooth_MKL COCG_Di_Smooth_MKL
#define PRECONDITIONER PRECONDITIONER_DI
#define COCG_LLT_SMOOTH_MKL_H_INCLUDED
#if !defined USE_MKL
#define NO_MKL_STUBS
#endif
#include "COCG_LLT_Smooth_MKL.cpp"

