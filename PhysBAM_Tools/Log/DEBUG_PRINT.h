//#####################################################################
// Copyright 2006, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Macro PHYSBAM_DEBUG_PRINT
//#####################################################################
#ifndef __DEBUG_PRINT__
#define __DEBUG_PRINT__

#ifdef WIN32

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#define PHYSBAM_DEBUG_PRINT(prefix,...) PHYSBAM_WARNING("Debug Print is not implemented on Windows.")

#else

namespace PhysBAM{

// see http://groups.google.com/group/comp.std.c/browse_thread/thread/7cdd9f2984c15e3e/effa7cda7c378dac%23effa7cda7c378dac for macro explanations

#define PHYSBAM_CAT(a,...) PHYSBAM_PRIMITIVE_CAT(a,__VA_ARGS__)
#define PHYSBAM_PRIMITIVE_CAT(a,...) a ## __VA_ARGS__
#define PHYSBAM_SPLIT(i,...) PHYSBAM_PRIMITIVE_CAT(PHYSBAM_SPLIT_,i)(__VA_ARGS__)
#define PHYSBAM_SPLIT_0(a,...) a
#define PHYSBAM_SPLIT_1(a,...) __VA_ARGS__
#define PHYSBAM_COMMA() ,
#define PHYSBAM_REM(...) __VA_ARGS__
#define PHYSBAM_SIZE(...) PHYSBAM_SPLIT(0,PHYSBAM_SPLIT(1,PHYSBAM_SIZE_A(PHYSBAM_COMMA,PHYSBAM_REM(__VA_ARGS__)),,))
#define PHYSBAM_SIZE_A(_,im) PHYSBAM_SIZE_B(im, _()15, _()14, _()13, _()12, _()11, _()10, _()9, _()8, _()7, _() 6, _() 5, _() 4, _() 3, _() 2, _() 1)
#define PHYSBAM_SIZE_B(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,_,...) _
#define PHYSBAM_OVERLOAD(prefix,...) PHYSBAM_CAT(prefix,PHYSBAM_SIZE(__VA_ARGS__))

#define PHYSBAM_MAP_1(f,a1)                                                  f(a1)
#define PHYSBAM_MAP_2(f,a1,a2)                                               f(a1),f(a2)
#define PHYSBAM_MAP_3(f,a1,a2,a3)                                            f(a1),f(a2),f(a3)
#define PHYSBAM_MAP_4(f,a1,a2,a3,a4)                                         f(a1),f(a2),f(a3),f(a4)
#define PHYSBAM_MAP_5(f,a1,a2,a3,a4,a5)                                      f(a1),f(a2),f(a3),f(a4),f(a5)
#define PHYSBAM_MAP_6(f,a1,a2,a3,a4,a5,a6)                                   f(a1),f(a2),f(a3),f(a4),f(a5),f(a6)
#define PHYSBAM_MAP_7(f,a1,a2,a3,a4,a5,a6,a7)                                f(a1),f(a2),f(a3),f(a4),f(a5),f(a6),f(a7)
#define PHYSBAM_MAP_8(f,a1,a2,a3,a4,a5,a6,a7,a8)                             f(a1),f(a2),f(a3),f(a4),f(a5),f(a6),f(a7),f(a8)
#define PHYSBAM_MAP_9(f,a1,a2,a3,a4,a5,a6,a7,a8,a9)                          f(a1),f(a2),f(a3),f(a4),f(a5),f(a6),f(a7),f(a8),f(a9)
#define PHYSBAM_MAP_10(f,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)                     f(a1),f(a2),f(a3),f(a4),f(a5),f(a6),f(a7),f(a8),f(a9),f(a10)
#define PHYSBAM_MAP_11(f,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11)                 f(a1),f(a2),f(a3),f(a4),f(a5),f(a6),f(a7),f(a8),f(a9),f(a10),f(a11)
#define PHYSBAM_MAP_12(f,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)             f(a1),f(a2),f(a3),f(a4),f(a5),f(a6),f(a7),f(a8),f(a9),f(a10),f(a11),f(a12)
#define PHYSBAM_MAP_13(f,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13)         f(a1),f(a2),f(a3),f(a4),f(a5),f(a6),f(a7),f(a8),f(a9),f(a10),f(a11),f(a12),f(a13)
#define PHYSBAM_MAP_14(f,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14)     f(a1),f(a2),f(a3),f(a4),f(a5),f(a6),f(a7),f(a8),f(a9),f(a10),f(a11),f(a12),f(a13),f(a14)
#define PHYSBAM_MAP_15(f,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15) f(a1),f(a2),f(a3),f(a4),f(a5),f(a6),f(a7),f(a8),f(a9),f(a10),f(a11),f(a12),f(a13),f(a14),f(a15)
#define PHYSBAM_MAP(f,...) PHYSBAM_OVERLOAD(PHYSBAM_MAP_,__VA_ARGS__)(f,__VA_ARGS__)

void Debug_Print_Helper(const char* prefix,...);

#define PHYSBAM_DEBUG_PRINT_HELPER(a) #a,PhysBAM::STRING_UTILITIES::Value_To_String(a).c_str()

#define PHYSBAM_DEBUG_PRINT(prefix,...) PhysBAM::Debug_Print_Helper(prefix,PHYSBAM_MAP(PHYSBAM_DEBUG_PRINT_HELPER,__VA_ARGS__),(char*)0)

}
#endif
#endif
