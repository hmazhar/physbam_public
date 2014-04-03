//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PHYSBAM_OVERRIDE__
#define __PHYSBAM_OVERRIDE__

#ifdef WIN32
#  define PHYSBAM_OVERRIDE override
#  define PHYSBAM_SEALED sealed
#  define PHYSBAM_UNUSED
#  define PHYSBAM_NORETURN(declaration) __declspec(noreturn) declaration
#  define PHYSBAM_ALWAYS_INLINE
#  define PHYSBAM_FLATTEN
#else
#  define PHYSBAM_OVERRIDE
#  define PHYSBAM_SEALED
#  define PHYSBAM_UNUSED __attribute__ ((unused))
#  define PHYSBAM_NORETURN(declaration) declaration __attribute__ ((noreturn))
#  ifdef NDEBUG
#    define PHYSBAM_ALWAYS_INLINE __attribute__ ((always_inline))
#    define PHYSBAM_FLATTEN __attribute__ ((flatten))
#  else
#    define PHYSBAM_ALWAYS_INLINE
#    define PHYSBAM_FLATTEN
#  endif
#endif

#endif
