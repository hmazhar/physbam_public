//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_GREEN_GRID_3D
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __RED_GREEN_GRID_3D__
#define __RED_GREEN_GRID_3D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Red_Green/RED_TETRAHEDRON.h>
namespace PhysBAM{

template<class T> class OCTREE_GRID;

template<class T>
class RED_GREEN_GRID_3D:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
public:
    typedef T SCALAR;
    typedef TV VECTOR_T;

    GRID<TV> uniform_grid; // grid points include nodes, edge midpoints, and tetrahedron centers
    ARRAY<RED_TETRAHEDRON<T>*,VECTOR<int,3> > elements; // only entries corresponding to tetrahedron centers are nonzero
    int number_of_cells;
    int number_of_nodes;
    int maximum_depth;
private:
    mutable bool node_locations_up_to_date;mutable ARRAY<VECTOR<T,3> > node_locations;
public:

    RED_GREEN_GRID_3D()
        :number_of_cells(0),number_of_nodes(0),maximum_depth(0)
    {
        Tree_Topology_Changed();
    }
    
    ~RED_GREEN_GRID_3D()
    {Clean_Memory();}
    
    void Clean_Memory()
    {number_of_cells=number_of_nodes=maximum_depth=0;
    for(int i=1;i<=uniform_grid.counts.x;i++)for(int j=1;j<=uniform_grid.counts.y;j++)for(int ij=1;ij<=uniform_grid.counts.z;ij++)if(elements(i,j,ij)){
        RED_CHILDREN_3D<T>* owner=elements(i,j,ij)->owner;owner->Clean_Memory();delete owner;}
    elements.Clean_Memory();}

    static RED_GREEN_GRID_3D<T>* Create()
    {return new RED_GREEN_GRID_3D<T>;}

    void Tree_Topology_Changed()
    {node_locations.Clean_Memory();node_locations_up_to_date=false;}

    void Set_Maximum_Depth(const int maximum_depth_input)
    {maximum_depth=maximum_depth_input;}

    ARRAY<VECTOR<T,3> >& Node_Locations() const
    {if(!node_locations_up_to_date){Calculate_Node_Locations(node_locations);node_locations_up_to_date=true;}
    return node_locations;}

    bool Outside(const VECTOR<T,3>& location) const
    {return uniform_grid.Outside(location);}

    template<class TV>
    TV Interpolate_Nodes(const ARRAY<TV>& u,const VECTOR<T,3>& X) const
    {Node_Locations();const RED_TETRAHEDRON<T>* red_leaf=Red_Leaf_Tetrahedron(X);int i,j,k,l;
        if(red_leaf->Has_Green_Children()) red_leaf->green_children->elements(red_leaf->green_children->Green_Leaf_Tetrahedron(X,node_locations)).Get(i,j,k,l);
    else{i=red_leaf->Node(0);j=red_leaf->Node(1);k=red_leaf->Node(2);l=red_leaf->Node(3);}
    VECTOR<T,3> weights=TETRAHEDRON<T>::Barycentric_Coordinates(X,node_locations(i),node_locations(j),node_locations(k),node_locations(l));
    return weights.x*u(i)+weights.y*u(j)+weights.z*u(k)+(1-weights.x-weights.y-weights.z)*u(l);}

    template<class TV>
    void Interpolate_Node_Values_To_New_Nodes(ARRAY<TV>& node_values,const int old_number_of_nodes)
    {for(int i=1;i<=uniform_grid.counts.x;i++)for(int j=1;j<=uniform_grid.counts.y;j++)for(int ij=1;ij<=uniform_grid.counts.z;ij++)if(elements(i,j,ij))
        elements(i,j,ij)->Interpolate_Node_Values_To_New_Nodes(node_values,old_number_of_nodes);}
    
    template<class TH>
    bool Refine_One_Level(TH* helper,bool (*refinement_criteria)(TH* helper,const RED_TETRAHEDRON<T>* tetrahedron))
    {bool refinement=false;int old_number_of_cells=number_of_cells;
    for(int i=1;i<=uniform_grid.counts.x;i++)for(int j=1;j<=uniform_grid.counts.y;j++)for(int ij=1;ij<=uniform_grid.counts.z;ij++)if(elements(i,j,ij))
        if(elements(i,j,ij)->Refine(number_of_cells,number_of_nodes,old_number_of_cells,*this,helper,refinement_criteria))refinement=true;
    for(int i=1;i<=uniform_grid.counts.x;i++)for(int j=1;j<=uniform_grid.counts.y;j++)for(int ij=1;ij<=uniform_grid.counts.z;ij++)if(elements(i,j,ij))
        elements(i,j,ij)->Finalize_Green_Refinement(number_of_cells);
    return refinement;}

//#####################################################################
    void Initialize(const GRID<TV>& uniform_grid_input,const int maximum_depth_input);
private:
    void Initialize_Tetrahedrons_Across_Primal_Face(const int axis,const VECTOR<int,3>& face,const ARRAY<int,VECTOR<int,3> >& nodes);
public:
    void Build_Mesh(TETRAHEDRON_MESH& tetrahedron_mesh,const ARRAY<T>* phi,ARRAY<int>& cell_to_tetrahedron_mapping,ARRAY<int>* node_to_particle_mapping) const;
    void Calculate_Node_Locations(ARRAY<VECTOR<T,3> >& node_locations) const;
    void Compact_Array_Indices(ARRAY<int>* cell_mapping_array,ARRAY<int>* node_mapping_array);
    const RED_TETRAHEDRON<T>* Red_Leaf_Tetrahedron(const VECTOR<T,3>& location) const;
    const RED_TETRAHEDRON<T>* Root_Level_Tetrahedron(const VECTOR<T,3>& location) const;
    void Create_Overlay_Grid(OCTREE_GRID<T>& grid,const int number_of_ghost_cells,const bool use_nodes,const bool use_faces) const;
//#####################################################################
};
}
#endif
#endif
