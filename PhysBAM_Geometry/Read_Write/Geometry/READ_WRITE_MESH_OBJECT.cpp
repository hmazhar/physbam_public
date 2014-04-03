//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_MESH_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_MESH_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_HEXAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

void Register_Read_Write_Mesh_Object()
{
#define DIMENSION_READ_WRITE_HELPER(T,RW) \
    Read_Write<STRUCTURE<VECTOR<T,1> >,RW>::Register_Read_Write<MESH_OBJECT<VECTOR<T,1>,POINT_SIMPLEX_MESH> >(); \
    Read_Write<STRUCTURE<VECTOR<T,1> >,RW>::Register_Read_Write<MESH_OBJECT<VECTOR<T,1>,SEGMENT_MESH> >(); \
    Read_Write<STRUCTURE<VECTOR<T,2> >,RW>::Register_Read_Write<MESH_OBJECT<VECTOR<T,2>,SEGMENT_MESH> >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<MESH_OBJECT<VECTOR<T,3>,SEGMENT_MESH> >(); \
    Read_Write<STRUCTURE<VECTOR<T,2> >,RW>::Register_Read_Write<MESH_OBJECT<VECTOR<T,2>,TRIANGLE_MESH> >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<MESH_OBJECT<VECTOR<T,3>,TRIANGLE_MESH> >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<MESH_OBJECT<VECTOR<T,3>,TETRAHEDRON_MESH> >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<MESH_OBJECT<VECTOR<T,3>,HEXAHEDRON_MESH> >();

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
template<class RW,class TV,class T_MESH> void Read_Write<MESH_OBJECT<TV,T_MESH>,RW>::
Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
{
    MESH_OBJECT<TV,T_MESH>& object=dynamic_cast<MESH_OBJECT<TV,T_MESH>&>(structure_object);
    int size;
    Read_Binary<RW>(input,object.mesh,size);
    object.particles.array_collection->Clean_Memory(); // strip everything away except for position
    object.particles.array_collection->Resize(size);
    Read_Binary_Array<RW>(input,object.particles.X.Get_Array_Pointer(),size);
    if(object.mesh.elements.m){
        int min_index=ARRAYS_COMPUTATIONS::Min(object.mesh.elements.Flattened()),max_index=ARRAYS_COMPUTATIONS::Max(object.mesh.elements.Flattened());
        if(min_index<1) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Invalid vertex index %d",min_index));
        if(max_index>object.particles.array_collection->Size()) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Read invalid vertex index %d (particles.array_collection->Size() = %d)",max_index,object.particles.array_collection->Size())); 
        object.Update_Number_Nodes();}
}
//#####################################################################
// Function Read_Structure_Helper
//#####################################################################
template<class RW,class TV,class T_MESH> void Read_Write<MESH_OBJECT<TV,T_MESH>,RW>::
Read_Structure_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
{
    MESH_OBJECT<TV,T_MESH>& object=dynamic_cast<MESH_OBJECT<TV,T_MESH>&>(structure_object);
    object.Clean_Memory();Read_Binary<RW>(input,object.mesh);
}
//#####################################################################
// Function Write_Helper
//#####################################################################
template<class RW,class TV,class T_MESH> void Read_Write<MESH_OBJECT<TV,T_MESH>,RW>::
Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
{
    const MESH_OBJECT<TV,T_MESH>& object=dynamic_cast<const MESH_OBJECT<TV,T_MESH>&>(structure_object);
    if(object.mesh.number_nodes!=object.particles.array_collection->Size()) PHYSBAM_FATAL_ERROR("number_nodes mismatch");
    Write_Binary<RW>(output,object.mesh,object.particles.X);
}
//#####################################################################
// Function Write_Structure_Helper
//#####################################################################
template<class RW,class TV,class T_MESH> void Read_Write<MESH_OBJECT<TV,T_MESH>,RW>::
Write_Structure_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
{
    const MESH_OBJECT<TV,T_MESH>& object=dynamic_cast<const MESH_OBJECT<TV,T_MESH>&>(structure_object);
    Write_Binary<RW>(output,object.mesh);
}
//#####################################################################
#define DIMENSION_HELPER(T,RW) \
template class Read_Write<MESH_OBJECT<VECTOR<T,1>,POINT_SIMPLEX_MESH>,RW>; \
template class Read_Write<MESH_OBJECT<VECTOR<T,1>,SEGMENT_MESH>,RW>; \
template class Read_Write<MESH_OBJECT<VECTOR<T,2>,SEGMENT_MESH>,RW>; \
template class Read_Write<MESH_OBJECT<VECTOR<T,3>,SEGMENT_MESH>,RW>; \
template class Read_Write<MESH_OBJECT<VECTOR<T,2>,TRIANGLE_MESH>,RW>; \
template class Read_Write<MESH_OBJECT<VECTOR<T,3>,TRIANGLE_MESH>,RW>; \
template class Read_Write<MESH_OBJECT<VECTOR<T,3>,TETRAHEDRON_MESH>,RW>;

DIMENSION_HELPER(float,float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
DIMENSION_HELPER(float,double);
DIMENSION_HELPER(double,float);
DIMENSION_HELPER(double,double);
#endif
}
#endif
