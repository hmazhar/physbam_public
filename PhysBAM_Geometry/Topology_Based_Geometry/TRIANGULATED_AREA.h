//#####################################################################
// Copyright 2003-2008, Christopher Allocco, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Igor Neverov, Duc Nguyen, Andrew Selle, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLULATED_AREA
//##################################################################### 
#ifndef __TRIANGULATED_AREA__
#define __TRIANGULATED_AREA__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SPATIAL_ACCELERATION_FORWARD.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class T_input>
class TRIANGULATED_AREA:public MESH_OBJECT<VECTOR<T_input,2>,TRIANGLE_MESH>
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;
    typedef TRIANGLE_HIERARCHY_2D<T>* T_HIERARCHY;
public:
    typedef MESH_OBJECT<TV,TRIANGLE_MESH> BASE;
    using BASE::mesh;using BASE::particles;using BASE::bounding_box;using BASE::Refresh_Auxiliary_Structures;using BASE::Rescale;

    SEGMENTED_CURVE_2D<T>* segmented_curve;
    TRIANGLE_HIERARCHY_2D<T>* hierarchy;
    ARRAY<VECTOR<T,2> >* triangle_area_fractions;
    ARRAY<T>* triangle_areas;
    ARRAY<T>* nodal_areas;

    TRIANGULATED_AREA(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input);
    virtual ~TRIANGULATED_AREA();

    void Compute_Nodal_And_Triangle_Areas()
    {Compute_Nodal_Areas(true);}

    SEGMENT_MESH& Get_Segment_Mesh()
    {return mesh.Get_Segment_Mesh();}

    SEGMENTED_CURVE_2D<T>& Get_Boundary_Object()
    {if(!segmented_curve) Initialize_Segmented_Curve();
    return *segmented_curve;}

    TRIANGLE_2D<T> Get_Element(const int aggregate_id) const
    {return TRIANGLE_2D<T>(particles.X.Subset(mesh.elements(aggregate_id)));}

    T Element_Size(const int triangle) const
    {return Area(triangle);}

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name()
    {return "SIMPLICIAL_OBJECT<T,VECTOR<T,2>,2>";}

    virtual std::string Extension() const PHYSBAM_OVERRIDE {return Static_Extension();}
    static std::string Static_Extension()
    {return "tri2d";}

    T Total_Size() const
    {return Total_Area(); }

    T Size(const int triangle) const
    {return Area(triangle);}

    T Signed_Size(const int triangle) const
    {return Signed_Area(triangle);}

//#####################################################################
    void Clean_Memory() PHYSBAM_OVERRIDE;
    void Initialize_Hierarchy(const bool update_boxes=true); // creates and updates the boxes as well
    void Initialize_Square_Mesh_And_Particles(const GRID<TV>& grid,const bool reverse_triangles=false);
    void Initialize_Circle_Mesh_And_Particles(const T outer_radius,const T inner_radius,const int num_radial,const int num_tangential);
    void Initialize_Herring_Bone_Mesh_And_Particles(const GRID<TV>& grid);
    void Initialize_Equilateral_Mesh_And_Particles(const GRID<TV>& grid);
    int Inside(const TV& location,const T thickness_over_two=0) const;
    void Check_Signed_Area_And_Make_Consistent(const int triangle,const bool verbose);
    void Check_Signed_Areas_And_Make_Consistent(const bool verbose=true);
    TV Centroid(const int triangle) const;
    void Rescale(const T scaling_factor) PHYSBAM_OVERRIDE;
    void Rescale(const T scaling_x,const T scaling_y);
    T Area(const int triangle) const;
    T Signed_Area(const int triangle) const;
    T Minimum_Area(int* index=0) const;
    T Minimum_Signed_Area(int* index=0) const;   
    T Total_Area() const;
    T Minimum_Altitude(int* index=0) const;
    T Minimum_Edge_Length(int* index=0) const;
    T Maximum_Edge_Length(int* index=0) const;
    void Initialize_Segmented_Curve();
    void Inverted_Triangles(ARRAY<int>& inverted_triangles) const;
    T Area_Incident_On_A_Particle(const int particle_index);
    int Split_Node(const int particle_index,const TV& normal);
    int Split_Connected_Component(const int node);
    void Discard_Triangles_Outside_Implicit_Curve(IMPLICIT_OBJECT<TV>& implicit_curve);
    void Initialize_Triangle_Area_Fractions_From_Voronoi_Regions();
    void Compute_Triangle_Areas();
    void Compute_Nodal_Areas(bool save_triangle_areas=false);
    int Triangle_In_Direction_Uninverted(const int node,const TV& direction) const;
    int Triangle_Walk_Uninverted(const int start_node,const TV& dX) const;
    bool Fix_Pair_For_Delaunay(const int triangle1,const int triangle2);
private:
    void Refresh_Auxiliary_Structures_Helper() PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

