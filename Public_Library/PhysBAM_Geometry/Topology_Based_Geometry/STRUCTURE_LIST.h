//#####################################################################
// Copyright 2003-2009, Eran Guendelman, Geoffrey Irving, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STRUCTURE_LIST
//#####################################################################
#ifndef __STRUCTURE_LIST__
#define __STRUCTURE_LIST__

#include <PhysBAM_Tools/Data_Structures/DYNAMIC_LIST.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
namespace PhysBAM{

template<class TV,class ID>
class STRUCTURE_LIST: public DYNAMIC_LIST<STRUCTURE<TV>,ID>
{
    typedef typename TV::SCALAR T;
    typedef DYNAMIC_LIST<STRUCTURE<TV>,ID> BASE;
public:
    mutable ARRAY<std::string> names;

    STRUCTURE_LIST()
        :BASE()
    {}
};
}
#endif
