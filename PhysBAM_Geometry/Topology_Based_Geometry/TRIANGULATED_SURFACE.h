//#####################################################################
// Copyright 2002-2008, Chris Allocco, Robert Bridson, Kevin Der, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Neil Molino, Igor Neverov, Andrew Selle, Jonathan Su, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGULATED_SURFACE
//##################################################################### 
#ifndef __TRIANGULATED_SURFACE__
#define __TRIANGULATED_SURFACE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SPATIAL_ACCELERATION_FORWARD.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class TV> class IMPLICIT_OBJECT;

template<class T_input>
class TRIANGULATED_SURFACE:public MESH_OBJECT<VECTOR<T_input,3>,TRIANGLE_MESH>
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
    typedef TRIANGLE_HIERARCHY<T> T_HIERARCHY;
public:
    typedef MESH_OBJECT<TV,TRIANGLE_MESH> BASE;
    using BASE::mesh;using BASE::particles;using BASE::Discard_Valence_Zero_Particles_And_Renumber;using BASE::bounding_box;using BASE::Refresh_Auxiliary_Structures;
    using BASE::Update_Bounding_Box;using BASE::Rescale;

    ARRAY<TRIANGLE_3D<T> >* triangle_list;
    ARRAY<T>* segment_lengths;
    TRIANGLE_HIERARCHY<T>* hierarchy;
    PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >* particle_hierarchy;
    bool use_vertex_normals,avoid_normal_interpolation_across_sharp_edges;
    T normal_variance_threshold;
    ARRAY<TV>* vertex_normals;    
    ARRAY<VECTOR<TV,3> >* face_vertex_normals;    
    
    TRIANGULATED_SURFACE(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input);
    virtual ~TRIANGULATED_SURFACE();

    void Use_Vertex_Normals() // averaged on the faces
    {use_vertex_normals=true;}
    
    void Use_Face_Normals() // i.e. flat normals
    {use_vertex_normals=false;}

    TRIANGLE_3D<T> Get_Element(const int aggregate_id) const
    {return TRIANGLE_3D<T>(particles.X.Subset(mesh.elements(aggregate_id)));}
    
    TV Face_Normal(const int index) const
    {if(triangle_list) return (*triangle_list)(index).normal;
    int i,j,k;mesh.elements(index).Get(i,j,k);return TRIANGLE_3D<T>::Normal(particles.X(i),particles.X(j),particles.X(k));}
    
    TV Centroid(const int triangle) const
    {int i,j,k;mesh.elements(triangle).Get(i,j,k);
    return (T)one_third*(particles.X(i)+particles.X(j)+particles.X(k));}

    T Area(const int triangle) const
    {int i,j,k;mesh.elements(triangle).Get(i,j,k);
    return TRIANGLE_3D<T>::Area(particles.X(i),particles.X(j),particles.X(k));}

    T Signed_Area(const int triangle) const
    {return Area(triangle);}

    void Rescale(const T scaling_factor) PHYSBAM_OVERRIDE
    {Rescale(scaling_factor,scaling_factor,scaling_factor);}

    SEGMENT_MESH& Get_Segment_Mesh()
    {return mesh.Get_Segment_Mesh();}

    T Element_Size(const int triangle) const
    {return Area(triangle);}

    TV Closest_Point_On_Boundary(const TV& location,const T max_depth=0,const T thickness_over_2=0,int* closest_triangle=0,T* distance=0) const // without max_depth, this is slow
    {return Surface(location,max_depth,thickness_over_2,closest_triangle,distance);}

    bool Inside_Any_Simplex(const TV& location,int& triangle_id,const T thickness_over_two=0) const
    {return Inside_Any_Triangle(location,triangle_id,thickness_over_two);}

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name()
    {return "SIMPLICIAL_OBJECT<T,VECTOR<T,3>,2>";}

    virtual std::string Extension() const PHYSBAM_OVERRIDE {return Static_Extension();}
    static std::string Static_Extension()
    {return "tri";}

    T Total_Size() const
    { return Total_Area(); }

    T Size(const int triangle) const
    {return Area(triangle);}

    T Signed_Size(const int triangle) const
    {return Signed_Area(triangle);}

    TV Normal(const int aggregate) const
    {return Face_Normal(aggregate);}

//#####################################################################
    void Clean_Memory() PHYSBAM_OVERRIDE;
    void Initialize_Hierarchy(const bool update_boxes=true,const int triangles_per_group=0); // creates and updates the boxes as well
    void Initialize_Particle_Hierarchy(const INDIRECT_ARRAY<ARRAY_VIEW<TV> >& particle_subset_input,const bool update_boxes=true,const int particles_per_group=10);
    void Rescale(const T scaling_x,const T scaling_y,const T scaling_z);
    void Update_Triangle_List(); // updates the triangles assuming the particle particles.Xs are already updated
    void Update_Triangle_List(ARRAY_VIEW<const TV> X);
    void Initialize_Torus_Mesh_And_Particles(const int m,const int n,const T major_radius,const T minor_radius); // TODO(jontg) ...
    void Initialize_Cylinder_Mesh_And_Particles(const int m,const int n,const T length,const T radius,const bool create_caps=true); // TODO(jontg) ...
    void Initialize_Segment_Lengths();
    void Update_Vertex_Normals();
    TV Normal(const TV& location,const int aggregate) const;
    bool Inside(const TV& location,const T thickness_over_two=0) const;
    bool Inside_Relative_To_Triangle(const TV& location,const int triangle_index_for_ray_test,const T thickness_over_two=0) const;
    bool Inside_Using_Ray_Test(RAY<TV>& ray,const T thickness_over_two=0) const;
    bool Outside(const TV& location,const T thickness_over_two=0) const;
    bool Boundary(const TV& location,const T thickness_over_two) const;
    bool Inside_Any_Triangle(const TV& location,int& triangle_id,const T thickness_over_two=0) const;
    TV Surface(const TV& location,const T max_depth=0,const T thickness_over_2=0,int* closest_triangle=0,T* distance=0) const; // without max_depth, this is slow
    TV Oriented_Surface(const TV& location,const TV &normal,const T max_depth=0,const T thickness_over_2=0,int* closest_triangle=0,T* distance=0) const;
    T Signed_Solid_Angle_Of_Triangle_Web(const TV& location,int web_root_node) const;
    bool Check_For_Self_Intersection(const T thickness_over_2=0,const bool update_bounding_boxes=true,ARRAY<VECTOR<int,2> >* intersecting_segment_triangle_pairs=0);
    bool Find_First_Segment_Triangle_Intersection(const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2,const int max_coarsening_attempts=5,
        const bool update_bounding_boxes=true);
    bool Segment_Triangle_Intersection(const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2=0,const bool update_bounding_boxes=true,
        ARRAY<VECTOR<int,2> >* intersecting_segment_triangle_pairs=0);
    void Get_Triangles_Near_Edges(ARRAY<ARRAY<int> >& triangles_near_edges,const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2=0,
        const bool update_bounding_boxes=true);
    TV Centroid_Of_Neighbors(const int node) const;
    T Calculate_Signed_Distance(const TV& location,T thickness=0) const;
    void Linearly_Subdivide();
    void Loop_Subdivide();
    void Root_Three_Subdivide();
    T Total_Area() const;
    T Minimum_Angle(int* index=0) const;
    T Maximum_Angle(int* index=0) const;
    T Average_Minimum_Angle() const;
    T Average_Maximum_Angle() const;
    T Minimum_Edge_Length(int* index=0) const;
    T Maximum_Edge_Length(int* index=0) const;
    T Average_Edge_Length() const;
    T Maximum_Aspect_Ratio(int* index=0) const;
    T Average_Aspect_Ratio() const;
    T Minimum_Area(int* index=0) const;
    T Minimum_Altitude(int* index=0) const;
    T Maximum_Magnitude_Phi(const IMPLICIT_OBJECT<TV>& implicit_surface,int* index=0);
    void Make_Orientations_Consistent_With_Implicit_Surface(const IMPLICIT_OBJECT<TV>& implicit_surface);
    void Close_Surface(const bool merge_coincident_vertices,const T merge_coincident_vertices_threshold,const bool fill_holes,const bool verbose=false);
    void Remove_Degenerate_Triangles(const T area_threshold=(T)1e-8);
    TRIANGULATED_SURFACE* Create_Compact_Copy() const;
private:
    void Refresh_Auxiliary_Structures_Helper() PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
