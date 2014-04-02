//#####################################################################
// Copyright 2009, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header DEBUG_CAST
//#####################################################################
#ifndef __DEBUG_CAST__
#define __DEBUG_CAST__

#ifdef NDEBUG
#define debug_cast static_cast
#else
#define debug_cast dynamic_cast
#endif

#endif
