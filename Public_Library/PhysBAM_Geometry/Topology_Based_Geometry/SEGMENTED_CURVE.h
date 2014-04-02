//#####################################################################
// Copyright 2003-2009, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENTED_CURVE
//##################################################################### 
#ifndef __SEGMENTED_CURVE__
#define __SEGMENTED_CURVE__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_1D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SPATIAL_ACCELERATION_FORWARD.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class SEGMENTED_CURVE:public MESH_OBJECT<TV,SEGMENT_MESH>
{
    typedef typename TV::SCALAR T;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::SEGMENT T_SEGMENT;
public:
    typedef MESH_OBJECT<TV,SEGMENT_MESH> BASE;typedef SEGMENT_HIERARCHY<TV> T_HIERARCHY;
    using BASE::mesh;using BASE::particles;using BASE::bounding_box;using BASE::Update_Bounding_Box;

    SEGMENT_HIERARCHY<TV>* hierarchy;
    ARRAY<T_SEGMENT>* segment_list;
    POINT_SIMPLICES_1D<T>* point_simplices_1d;

    SEGMENTED_CURVE(SEGMENT_MESH& segment_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input);
    virtual ~SEGMENTED_CURVE();

    T_SEGMENT Get_Element(const int aggregate_id) const
    {return T_SEGMENT(particles.X(mesh.elements(aggregate_id)(1)),particles.X(mesh.elements(aggregate_id)(2)));}

    void Rescale(const T scaling_factor) PHYSBAM_OVERRIDE
    {Rescale(scaling_factor*TV::All_Ones_Vector());}

    SEGMENT_MESH& Get_Segment_Mesh()
    {return mesh;}

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name()
    {return STRING_UTILITIES::string_sprintf("SIMPLICIAL_OBJECT<T,VECTOR<T,%d>,1>",TV::dimension);}

    virtual std::string Extension() const PHYSBAM_OVERRIDE {return Static_Extension();}
    static std::string Static_Extension()
    {return TV::dimension==2?"curve2d":(TV::dimension==1?"curve1d":"curve");}

    TV Closest_Point_On_Boundary(const TV& location,const T max_depth=0,const T thickness_over_2=0,int* closest_segment=0,T* distance=0) const
    {return Closest_Point_On_Curve(location,thickness_over_2,closest_segment,distance);} // max_depth currently unused

    T Total_Size() const
    {return Total_Length(); }

    T Size(const int segment) const
    {int node1=mesh.elements(segment)(1),node2=mesh.elements(segment)(2);
    return (particles.X(node1)-particles.X(node2)).Magnitude();}

    T Signed_Size(const int segment) const
    {return Size(segment);}

//#####################################################################
    POINT_SIMPLICES_1D<T>& Get_Boundary_Object();
    void Clean_Memory() PHYSBAM_OVERRIDE;
    void Update_Segment_List();
    void Initialize_Hierarchy(const bool update_boxes=true); // creates and updates the boxes as well
    void Initialize_Straight_Mesh_And_Particles(const GRID<VECTOR<T,1> >& grid);
    void Initialize_Circle_Mesh_And_Particles(const int m,const T radius);
    TV Closest_Point_On_Curve(const TV& location,T thickness_over_two=0,int* closest_segment=0,T* distance=0) const;
    void Rescale(const TV& scaling);
    T Average_Edge_Length() const;
    T Total_Length() const;
private:
    void Refresh_Auxiliary_Structures_Helper() PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
