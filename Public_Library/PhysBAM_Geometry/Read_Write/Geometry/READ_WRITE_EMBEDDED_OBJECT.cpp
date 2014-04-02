//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_EMBEDDED_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD_SUBSET.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_EMBEDDED_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
namespace PhysBAM{

void Register_Read_Write_Embedded_Object()
{
#define DIMENSION_READ_WRITE_HELPER(T,RW) \
    Read_Write<STRUCTURE<VECTOR<T,2> >,RW>::Register_Read_Write<EMBEDDED_OBJECT<VECTOR<T,2>,2> >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<EMBEDDED_OBJECT<VECTOR<T,3>,2> >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<EMBEDDED_OBJECT<VECTOR<T,3>,3> >(); \

#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
#define RW_HELPER(RW) \
    DIMENSION_READ_WRITE_HELPER(float,RW)
#else
#define RW_HELPER(RW) \
    DIMENSION_READ_WRITE_HELPER(float,RW) \
    DIMENSION_READ_WRITE_HELPER(double,RW)
#endif

    RW_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    RW_HELPER(double);
#endif
}

//#####################################################################
// Function Read_Helper
//#####################################################################
template<class RW,class TV,int d> void Read_Write<EMBEDDED_OBJECT<TV,d>,RW>::
Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
{
    EMBEDDED_OBJECT<TV,d>& object=dynamic_cast<EMBEDDED_OBJECT<TV,d>&>(structure_object);
    object.Clean_Memory(); // reads entire simplicial_object instead of just simplicial_object.mesh
    int backward_compatible;
    Read_Binary<RW>(input,object.embedded_particles,backward_compatible,object.parent_particles,object.interpolation_fraction,object.simplicial_object,object.embedded_object.mesh,backward_compatible);
    Read_Binary<RW>(input,object.node_in_simplex_is_material,object.interpolation_fraction_threshold,object.orientation_index);
}
//#####################################################################
// Function Read_Structure_Helper
//#####################################################################
template<class RW,class TV,int d> void Read_Write<EMBEDDED_OBJECT<TV,d>,RW>::
Read_Structure_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
{
    EMBEDDED_OBJECT<TV,d>& object=dynamic_cast<EMBEDDED_OBJECT<TV,d>&>(structure_object);
    object.Clean_Memory();
    int backward_compatible;
    Read_Binary<RW>(input,object.embedded_particles,backward_compatible,object.parent_particles,object.interpolation_fraction,object.simplicial_object.mesh,object.embedded_object.mesh,backward_compatible);
    Read_Binary<RW>(input,object.node_in_simplex_is_material,object.interpolation_fraction_threshold,object.orientation_index);
}
//#####################################################################
// Function Write_Helper
//#####################################################################
template<class RW,class TV,int d> void Read_Write<EMBEDDED_OBJECT<TV,d>,RW>::
Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
{
    const EMBEDDED_OBJECT<TV,d>& object=dynamic_cast<const EMBEDDED_OBJECT<TV,d>&>(structure_object);
    Write_Binary<RW>(output,object.embedded_particles,2,object.parent_particles,object.interpolation_fraction,object.simplicial_object,object.embedded_object.mesh,d+1,object.node_in_simplex_is_material);
    Write_Binary<RW>(output,object.interpolation_fraction_threshold,object.orientation_index);
}
//#####################################################################
// Function Write_Structure_Helper
//#####################################################################
template<class RW,class TV,int d> void Read_Write<EMBEDDED_OBJECT<TV,d>,RW>::
Write_Structure_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
{
    const EMBEDDED_OBJECT<TV,d>& object=dynamic_cast<const EMBEDDED_OBJECT<TV,d>&>(structure_object);
    Write_Binary<RW>(output,object.embedded_particles,2,object.parent_particles,object.interpolation_fraction,object.simplicial_object.mesh,object.embedded_object.mesh,d+1,object.node_in_simplex_is_material);
    Write_Binary<RW>(output,object.interpolation_fraction_threshold,object.orientation_index);
}

#define DIMENSION_HELPER(T,RW) \
    template class Read_Write<EMBEDDED_OBJECT<VECTOR<T,2>,2>,RW>; \
    template class Read_Write<EMBEDDED_OBJECT<VECTOR<T,3>,2>,RW>; \
    template class Read_Write<EMBEDDED_OBJECT<VECTOR<T,3>,3>,RW>;
//#####################################################################
DIMENSION_HELPER(float,float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
DIMENSION_HELPER(float,double);
DIMENSION_HELPER(double,float);
DIMENSION_HELPER(double,double);
#endif
}
#endif
