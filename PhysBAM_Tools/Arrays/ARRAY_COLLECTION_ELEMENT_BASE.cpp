//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION_ELEMENT_BASE.h>
using namespace PhysBAM;
ARRAY_COLLECTION_ELEMENT_BASE::
ARRAY_COLLECTION_ELEMENT_BASE()
    :id(ATTRIBUTE_ID()),owns_data(true)
{
}
ARRAY_COLLECTION_ELEMENT_BASE::
~ARRAY_COLLECTION_ELEMENT_BASE()
{
}
