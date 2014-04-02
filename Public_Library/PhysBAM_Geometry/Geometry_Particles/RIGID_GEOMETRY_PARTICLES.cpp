//#####################################################################
// Copyright 2006-2009, Michael Lentine, Craig Schroeder, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/RIGID_GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
namespace PhysBAM{
//#####################################################################
// Resize
//#####################################################################
template<class TV> void RIGID_GEOMETRY_PARTICLES<TV>::
Resize(const int new_size)
{
    for(int p=new_size+1;p<=array_collection->Size();p++)
        if(rigid_geometry(p)) Remove_Geometry(p);
    array_collection->Resize(new_size);
}
//#####################################################################
// Remove_Geometry
//#####################################################################
template<class TV> void RIGID_GEOMETRY_PARTICLES<TV>::
Remove_Geometry(const int p)
{
    delete rigid_geometry(p);
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void RIGID_GEOMETRY_PARTICLES<TV>::
Clean_Memory()
{
    for(int p=1;p<=array_collection->Size();p++) if(rigid_geometry(p)) Remove_Geometry(p);
    array_collection->Clean_Memory();
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void RIGID_GEOMETRY_PARTICLES<TV>::
Delete_All_Particles()
{
    for(int p=1;p<=array_collection->Size();p++) if(rigid_geometry(p)) Remove_Geometry(p);
    array_collection->Delete_All_Elements();
}
template class RIGID_GEOMETRY_PARTICLES<VECTOR<float,1> >;
template class RIGID_GEOMETRY_PARTICLES<VECTOR<float,2> >;
template class RIGID_GEOMETRY_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_GEOMETRY_PARTICLES<VECTOR<double,1> >;
template class RIGID_GEOMETRY_PARTICLES<VECTOR<double,2> >;
template class RIGID_GEOMETRY_PARTICLES<VECTOR<double,3> >;
#endif
}
