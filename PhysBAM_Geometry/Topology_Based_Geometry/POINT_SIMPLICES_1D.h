//#####################################################################
// Copyright 2007-2009, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINT_SIMPLICES_1D
//##################################################################### 
#ifndef __POINT_SIMPLICES_1D__
#define __POINT_SIMPLICES_1D__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SPATIAL_ACCELERATION_FORWARD.h>
#include <PhysBAM_Geometry/Topology/POINT_SIMPLEX_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
namespace PhysBAM{

template<class T_input>
class POINT_SIMPLICES_1D:public MESH_OBJECT<VECTOR<T_input,1>,POINT_SIMPLEX_MESH>
{
    typedef T_input T;
    STATIC_ASSERT(IS_FLOAT_OR_DOUBLE<T>::value);
    typedef VECTOR<T,1> TV;
public:
    typedef MESH_OBJECT<TV,POINT_SIMPLEX_MESH> BASE;
    using BASE::mesh;using BASE::particles;using BASE::bounding_box;

    ARRAY<int> indices;
    ARRAY<POINT_SIMPLEX_1D<T> >* point_simplex_list;
    PARTICLE_PARTITION<TV>* particle_partition;
    PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >* hierarchy;
    int number_point_simplices;

    POINT_SIMPLICES_1D(POINT_SIMPLEX_MESH& point_simplex_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input);

    void Clean_Memory() PHYSBAM_OVERRIDE
    {MESH_OBJECT<TV,POINT_SIMPLEX_MESH>::Clean_Memory();}

    bool Boundary(const TV& location,T thickness_over_two=0) const
    {PHYSBAM_NOT_IMPLEMENTED();}

    TV Normal(const int aggregate) const
    {PHYSBAM_NOT_IMPLEMENTED();}

    TV Normal(const TV& location,const int aggregate) const
    {PHYSBAM_NOT_IMPLEMENTED();}

    POINT_SIMPLEX_1D<T> Get_Element(const int aggregate_id) const
    {int node1;mesh.elements(aggregate_id).Get(node1);
    return POINT_SIMPLEX_1D<T>(particles.X(node1),mesh.directions(node1));}

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name()
    {return STRING_UTILITIES::string_sprintf("SIMPLICIAL_OBJECT<T,VECTOR<T,%d>,0>",TV::dimension);}

    virtual std::string Extension() const PHYSBAM_OVERRIDE {return Static_Extension();}
    static std::string Static_Extension() 
    {return "ptsimp";}

    void Update_Point_Simplex_List()
    {Update_Point_Simplex_List(particles.X);}

    template<class T_ARRAY_TV> void Update_Point_Simplex_List(const T_ARRAY_TV& X)
    {STATIC_ASSERT((IS_SAME<TV,typename T_ARRAY_TV::ELEMENT>::value));
    int number_point_simplices=mesh.elements.m;
    if(!point_simplex_list) point_simplex_list=new ARRAY<POINT_SIMPLEX_1D<T> >(number_point_simplices);else point_simplex_list->Resize(number_point_simplices);
    for(int k=1;k<=number_point_simplices;k++){
        int node1;mesh.elements(k).Get(node1);
        (*point_simplex_list)(k).x1=X(node1);
        (*point_simplex_list)(k).direction=mesh.directions(node1);}}

    SEGMENT_MESH& Get_Segment_Mesh()
    {return mesh.Get_Segment_Mesh();}

//#####################################################################
    TV Closest_Point_On_Boundary(const TV& location,const T max_depth=0,const T thickness_over_two=0,int* closest_point_simplex=0,T* distance=0) const;
    T Calculate_Signed_Distance(const TV& location,T thickness=0) const;
    bool Inside(const TV& location,T thickness_over_two=0) const;
    bool Outside(const TV& location,T thickness_over_two=0) const;
    bool Inside_Any_Simplex(const TV& point,int& point_simplex_id,const T thickness_over_two=0) const;
    void Initialize_Hierarchy(const bool update_boxes=true,const int elements_per_group=0)
    {
        indices.Resize(particles.X.Size()); for(int i=1;i<=particles.X.Size();++i) indices(i)=i;
        delete hierarchy; hierarchy=new PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >(INDIRECT_ARRAY<ARRAY_VIEW<TV> >(particles.X,indices),update_boxes,elements_per_group);
    }
    void Initialize_Particle_Partition(const VECTOR<int,1>& counts){PHYSBAM_NOT_IMPLEMENTED();}
private:
    void Refresh_Auxiliary_Structures_Helper() PHYSBAM_OVERRIDE {if(point_simplex_list) Update_Point_Simplex_List();}
//#####################################################################
};
}
#endif
