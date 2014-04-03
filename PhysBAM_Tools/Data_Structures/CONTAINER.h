//#####################################################################
// Copyright 2007, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONTAINER
//#####################################################################
#ifndef __CONTAINER__
#define __CONTAINER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
namespace PhysBAM{

template<class ID1>
class CONTAINER
{
public:
    typedef ID1 ID;
    ID id;
};

template<class ID1,class ID2>
class TOP_LEVEL_CONTAINER
{
public:
    typedef ID1 ID;
    ID id;
    ARRAY<ID2> containers;
};
}
#endif
