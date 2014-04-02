//#####################################################################
// Copyright 2004, Ron Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_OBJECT_ACCELERATION_PRIMITIVE
//#####################################################################
#ifndef __RENDERING_OBJECT_ACCELERATION_PRIMITIVE__
#define __RENDERING_OBJECT_ACCELERATION_PRIMITIVE__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
namespace PhysBAM{

template<class T> class RENDERING_OBJECT;

template<class T>
class RENDERING_OBJECT_ACCELERATION_PRIMITIVE
{
public:
    RANGE<VECTOR<T,3> > world_bounding_box;
    const RENDERING_OBJECT<T>* object;
    int aggregate_id;
    int operation;
    T hit_t;
    int hit_aggregate_id;

    RENDERING_OBJECT_ACCELERATION_PRIMITIVE(const RANGE<VECTOR<T,3> >& bounding_box_in,const RENDERING_OBJECT<T>* object_in,const int aggregate_id_in)
        :world_bounding_box(bounding_box_in),object(object_in),aggregate_id(aggregate_id_in),operation(0),hit_aggregate_id(0)
    {}

    RENDERING_OBJECT_ACCELERATION_PRIMITIVE()
    {}

//#####################################################################
};
}
#endif
