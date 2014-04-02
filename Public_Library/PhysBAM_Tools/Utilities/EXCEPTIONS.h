//#####################################################################
// Copyright 2008, Geoffrey Irving, Ryan Kautzman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// File EXCEPTIONS
//#####################################################################
#ifndef __EXCEPTIONS__
#define __EXCEPTIONS__

#include <stdexcept>
namespace PhysBAM{

// Note: destructors must be in .cpp to avoid shared library name lookup issues

#define PHYSBAM_SIMPLE_EXCEPTION(ERROR,BASE) \
    struct ERROR:public BASE                 \
    {                                        \
        ERROR(const std::string& message)    \
            :BASE(message)                   \
        {}                                   \
                                             \
        virtual ~ERROR() throw ();           \
    };

PHYSBAM_SIMPLE_EXCEPTION(PHYSBAM_ERROR,std::runtime_error) // base class for all physbam exceptions

PHYSBAM_SIMPLE_EXCEPTION(READ_ERROR,PHYSBAM_ERROR)
PHYSBAM_SIMPLE_EXCEPTION(FILESYSTEM_ERROR,PHYSBAM_ERROR)
PHYSBAM_SIMPLE_EXCEPTION(LOOKUP_ERROR,PHYSBAM_ERROR)
    PHYSBAM_SIMPLE_EXCEPTION(INDEX_ERROR,LOOKUP_ERROR)
    PHYSBAM_SIMPLE_EXCEPTION(KEY_ERROR,LOOKUP_ERROR)
PHYSBAM_SIMPLE_EXCEPTION(TYPE_ERROR,PHYSBAM_ERROR)
PHYSBAM_SIMPLE_EXCEPTION(VALUE_ERROR,PHYSBAM_ERROR)
PHYSBAM_SIMPLE_EXCEPTION(NOT_IMPLEMENTED_ERROR,PHYSBAM_ERROR)
PHYSBAM_SIMPLE_EXCEPTION(ASSERTION_ERROR,PHYSBAM_ERROR)
PHYSBAM_SIMPLE_EXCEPTION(FLOATING_POINT_ERROR,PHYSBAM_ERROR)

#undef PHYSBAM_SIMPLE_EXCEPTION

}
#endif
