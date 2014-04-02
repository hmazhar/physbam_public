//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Robert Bridson, Ron Fedkiw, Eilene Hao, Geoffrey Irving, Neil Molino, Avi Robinson-Mosher, Mike Rodgers, Andrew Selle, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRALIZED_VOLUME
//#####################################################################
#ifndef __TETRAHEDRALIZED_VOLUME__
#define __TETRAHEDRALIZED_VOLUME__

#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
namespace PhysBAM{

template<class TV> class IMPLICIT_OBJECT;
template<class T> class TETRAHEDRON_HIERARCHY;
template<class TV> class GRID;

template<class T_input>
class TETRAHEDRALIZED_VOLUME:public MESH_OBJECT<VECTOR<T_input,3>,TETRAHEDRON_MESH>
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef MESH_OBJECT<TV,TETRAHEDRON_MESH> BASE;
    using BASE::mesh;using BASE::particles;using BASE::bounding_box;using BASE::Update_Bounding_Box;

    ARRAY<TETRAHEDRON<T> >* tetrahedron_list;
    TRIANGULATED_SURFACE<T>* triangulated_surface;
    TETRAHEDRON_HIERARCHY<T>* hierarchy;
    ARRAY<T>* tetrahedron_volumes;
    ARRAY<T>* nodal_volumes;

    TETRAHEDRALIZED_VOLUME(TETRAHEDRON_MESH& tetrahedron_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input);
    virtual ~TETRAHEDRALIZED_VOLUME();

    void Compute_Nodal_And_Tetrahedron_Volumes()
    {Compute_Nodal_Volumes(true);}

    SEGMENT_MESH& Get_Segment_Mesh()
    {return mesh.Get_Segment_Mesh();}

    TRIANGULATED_SURFACE<T>& Get_Boundary_Object()
    {if(!triangulated_surface) Initialize_Triangulated_Surface();
    return *triangulated_surface;}

    TETRAHEDRON<T> Get_Element(const int aggregate_id) const
    {return TETRAHEDRON<T>(particles.X.Subset(mesh.elements(aggregate_id)));}

    T Element_Size(const int tetrahedron) const
    {return Volume(tetrahedron);}

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name()
    {return "SIMPLICIAL_OBJECT<T,VECTOR<T,3>,3>";}

    virtual std::string Extension() const PHYSBAM_OVERRIDE {return Static_Extension();}
    static std::string Static_Extension()
    {return "tet";}

    T Total_Size() const
    {return Total_Volume(); }

    T Size(const int tetrahedron) const
    {return Volume(tetrahedron);}

    T Signed_Size(const int tetrahedron) const
    {return Signed_Volume(tetrahedron);}

//#####################################################################
    void Clean_Memory() PHYSBAM_OVERRIDE;
    void Update_Tetrahedron_List(); // updates the tets assuming the particle positions are already updated
    void Initialize_Hierarchy(const bool update_boxes=true); // creates and updates the boxes as well
    void Initialize_Octahedron_Mesh_And_Particles(const GRID<TV>& grid);
    void Initialize_Cube_Mesh_And_Particles(const GRID<TV>& grid);
    void Initialize_Prismatic_Cube_Mesh_And_Particles(const GRID<TV>& grid);
    void Check_Signed_Volumes_And_Make_Consistent(bool verbose=true);
    void Initialize_Triangulated_Surface();
    T Minimum_Volume(int* index=0) const;
    T Minimum_Signed_Volume(int* index=0) const;
    T Total_Volume() const;
    T Minimum_Angle(int* index=0) const;
    T Maximum_Angle(int* index=0) const;
    T Minimum_Altitude(int* index=0) const;
    T Minimum_Edge_Length(int* index=0) const;
    T Maximum_Edge_Length(int* index=0) const;
    T Maximum_Aspect_Ratio(int* index=0) const;
    T Maximum_Interior_Aspect_Ratio(int* index=0);
    T Maximum_Boundary_Aspect_Ratio(int* index=0);
    T Average_Aspect_Ratio();
    T Average_Interior_Aspect_Ratio();
    T Average_Boundary_Aspect_Ratio();
    T Minimum_Dihedral_Angle(int* index=0) const;
    T Maximum_Dihedral_Angle(int* index=0) const;
    T Maximum_Edge_Length(int* index=0);
    T Minimum_Edge_Length(int* index=0);
    void Advance_Interior_Laplacian_Smoothing();
    TV Centroid(const int tetrahedron) const;
    void Rescale(const T scaling_factor) PHYSBAM_OVERRIDE;
    void Rescale(const T scaling_x,const T scaling_y,const T scaling_z);
    T Volume(const int tetrahedron) const;
    T Signed_Volume(const int tetrahedron) const;
    TV Centroid_Of_Neighbors(const int node) const;
    bool Completely_Inside_Box(const int tetrahedron,const RANGE<TV>& box) const;
    void Discard_Spikes_From_Adjacent_Elements(ARRAY<int>* deletion_list=0);
    void Discard_Spikes(ARRAY<int>* deletion_list=0);
    void Interior_Edges_With_Boundary_Nodes(ARRAY<VECTOR<int,2> >* deletion_list);
    void Inverted_Tetrahedrons(ARRAY<int>& inverted_tetrahedrons) const;
    int Inside(const TV& location,const T thickness_over_two=0) const;
    void Discard_Tetrahedrons_Outside_Implicit_Surface(IMPLICIT_OBJECT<TV>& implicit_surface);
    void Discard_Tetrahedrons_Outside_Implicit_Surface_Aggressive(IMPLICIT_OBJECT<TV>& implicit_surface);
    void Discard_Tetrahedrons_Outside_Implicit_Surface_Aggressive(IMPLICIT_OBJECT<TV>& implicit_surface,const ARRAY<RANGE<TV> >& bounding_boxes);
    T Maximum_Magnitude_Phi_On_Boundary(const IMPLICIT_OBJECT<TV>& implicit_surface,int* index=0);
    T Volume_Incident_On_A_Particle(const int particle_index);
    void Split_Along_Fracture_Plane(const PLANE<T>& plane,ARRAY<int>& particle_replicated);
    int Split_Node(const int particle_index,const TV& normal);
    int Split_Connected_Component(const int node);
    void Compute_Tetrahedron_Volumes();
    void Compute_Nodal_Volumes(bool save_tetrahedron_volumes=false);
private:
    void Refresh_Auxiliary_Structures_Helper() PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
