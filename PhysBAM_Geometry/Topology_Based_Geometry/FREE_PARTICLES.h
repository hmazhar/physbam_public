//#####################################################################
// Copyright 2006, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//######################################################################
// Class FREE_PARTICLES
//######################################################################
#ifndef __FREE_PARTICLES__
#define __FREE_PARTICLES__

#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
namespace PhysBAM{

template<class TV>
class FREE_PARTICLES:public STRUCTURE<TV>
{
    typedef typename TV::SCALAR T;
public:
    ARRAY<int> nodes;

    FREE_PARTICLES();
    ~FREE_PARTICLES();

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name()
    {return STRING_UTILITIES::string_sprintf("FREE_PARTICLES<T,VECTOR<T,%d> >",TV::m);}

    static FREE_PARTICLES* Create(GEOMETRY_PARTICLES<TV>& particles)
    {return Create();}

//######################################################################
    static FREE_PARTICLES* Create();
    void Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const PHYSBAM_OVERRIDE;
//######################################################################
};
}
#endif
