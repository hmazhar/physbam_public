//#####################################################################
// Copyright 2008, Geoffrey Irving, Ryan Kautzman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// File EXCEPTIONS
//#####################################################################
#include <PhysBAM_Tools/Utilities/EXCEPTIONS.h>
namespace PhysBAM{

#define INSTANTIATE(ERROR) \
    ERROR::~ERROR() throw () {}

INSTANTIATE(PHYSBAM_ERROR)

INSTANTIATE(READ_ERROR)
INSTANTIATE(FILESYSTEM_ERROR)
INSTANTIATE(LOOKUP_ERROR)
INSTANTIATE(INDEX_ERROR)
INSTANTIATE(KEY_ERROR)
INSTANTIATE(TYPE_ERROR)
INSTANTIATE(VALUE_ERROR)
INSTANTIATE(NOT_IMPLEMENTED_ERROR)
INSTANTIATE(ASSERTION_ERROR)
INSTANTIATE(FLOATING_POINT_ERROR)

}
