//#####################################################################
// Copyright 2006-2009, Michael Lentine, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GEOMETRY_PARTICLES
//#####################################################################
#ifndef __RIGID_GEOMETRY_PARTICLES__
#define __RIGID_GEOMETRY_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>

namespace PhysBAM{

template<class TV> class RIGID_GEOMETRY;
template<class TV> class STRUCTURE;
template<class TV> class COLLISION_GEOMETRY_COLLECTION;
template<class TV,class ID> class STRUCTURE_LIST;

template<class TV>
class RIGID_GEOMETRY_PARTICLES:public CLONEABLE<RIGID_GEOMETRY_PARTICLES<TV>,POINT_CLOUD<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<RIGID_GEOMETRY_PARTICLES<TV>,POINT_CLOUD<TV> > BASE;
public:
    using BASE::array_collection;using BASE::X;

    ARRAY_VIEW<RIGID_GEOMETRY<TV>*> rigid_geometry;
    ARRAY_VIEW<ROTATION<TV> > rotation;
    ARRAY_VIEW<TWIST<TV> > twist;
    ARRAY_VIEW<VECTOR<int,3> > structure_ids;

    RIGID_GEOMETRY_PARTICLES(ARRAY_COLLECTION* array_collection_input)
        :rigid_geometry(0,0),rotation(0,0),twist(0,0),structure_ids(0,0)
    {delete array_collection;array_collection=array_collection_input;Initialize_Array_Collection();}

    RIGID_GEOMETRY_PARTICLES()
        :rigid_geometry(0,0),rotation(0,0),twist(0,0),structure_ids(0,0)
    {Initialize_Array_Collection();}

    void Initialize_Array_Collection()
    {POINT_CLOUD<TV>::Initialize_Array_Collection();
    array_collection->Add_Array(ATTRIBUTE_ID_RIGID_GEOMETRY,&rigid_geometry);array_collection->Add_Array(ATTRIBUTE_ID_ROTATION,&rotation);
    array_collection->Add_Array(ATTRIBUTE_ID_TWIST,&twist);array_collection->Add_Array(ATTRIBUTE_ID_STRUCTURE_IDS,&structure_ids);}
    
    ~RIGID_GEOMETRY_PARTICLES()
    {}

//#####################################################################
    void Resize(const int new_size);
    void Remove_Geometry(const int p);
    void Clean_Memory();
    void Delete_All_Particles();
//#####################################################################
};
}
#endif
