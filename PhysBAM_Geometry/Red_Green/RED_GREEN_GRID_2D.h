//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_GREEN_GRID_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __RED_GREEN_GRID_2D__
#define __RED_GREEN_GRID_2D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Red_Green/RED_TRIANGLE.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
namespace PhysBAM{

template<class T>
class RED_GREEN_GRID_2D:public NONCOPYABLE
{
    typedef VECTOR<T,2> TV;
public:
    typedef T SCALAR;
    typedef TV VECTOR_T;

    GRID<TV> uniform_grid; // grid points include nodes, edge midpoints, and triangle centers
    ARRAY<RED_TRIANGLE<T>*,VECTOR<int,2> > elements; // only entries corresponding to triangle centers are nonzero
    int number_of_cells;
    int number_of_nodes;
    int maximum_depth;
    T max_edge_of_min_triangle;
private:
    mutable bool cell_pointer_from_index_up_to_date;mutable ARRAY<RED_TRIANGLE<T>*> cell_pointer_from_index;
    mutable bool node_locations_up_to_date;mutable ARRAY<VECTOR<T,2> > node_locations;
    mutable bool node_neighbors_up_to_date;mutable ARRAY<ARRAY<int> > node_neighbors;
    mutable bool triangle_mesh_up_to_date;mutable TRIANGLE_MESH triangle_mesh;mutable ARRAY<int> cell_to_triangle_mapping;
public:

    RED_GREEN_GRID_2D()
        :number_of_cells(0),number_of_nodes(0),maximum_depth(0)
    {
        Tree_Topology_Changed();
    }
    
    ~RED_GREEN_GRID_2D()
    {Clean_Memory();}
    
    void Clean_Memory()
    {number_of_cells=number_of_nodes=maximum_depth=0;
    for(int i=1;i<=uniform_grid.counts.x;i++)for(int j=1;j<=uniform_grid.counts.y;j++)if(elements(i,j)){
        RED_CHILDREN_2D<T>* owner=elements(i,j)->owner;owner->Clean_Memory();delete owner;}
    elements.Clean_Memory();}

    static RED_GREEN_GRID_2D<T>* Create()
    {return new RED_GREEN_GRID_2D<T>;}

    void Tree_Topology_Changed()
    {node_locations.Clean_Memory();node_locations_up_to_date=false;
    node_neighbors.Clean_Memory();node_neighbors_up_to_date=false;
    cell_pointer_from_index.Clean_Memory();cell_pointer_from_index_up_to_date=false;
    triangle_mesh.Clean_Memory();cell_to_triangle_mapping.Clean_Memory();triangle_mesh_up_to_date=false;}

    T Compute_Max_Edge_Of_Min_Triangle(const int proposed_maximum_depth) const
    {return uniform_grid.dX.x/(1<<(proposed_maximum_depth-2));}

    void Set_Maximum_Depth(const int maximum_depth_input)
    {maximum_depth=maximum_depth_input;max_edge_of_min_triangle=Compute_Max_Edge_Of_Min_Triangle(maximum_depth);}

    ARRAY<RED_TRIANGLE<T>*>& Cell_Pointer_From_Index()
    {if(!cell_pointer_from_index_up_to_date){Calculate_Cell_Pointer_From_Index(cell_pointer_from_index);cell_pointer_from_index_up_to_date=true;}
    return cell_pointer_from_index;}
    
    ARRAY<VECTOR<T,2> >& Node_Locations() const
    {if(!node_locations_up_to_date){Calculate_Node_Locations(node_locations);node_locations_up_to_date=true;}
    return node_locations;}

    ARRAY<ARRAY<int> >& Node_Neighbors() const
    {if(!node_neighbors_up_to_date){Calculate_Node_Neighbors(node_neighbors);node_neighbors_up_to_date=true;}
    return node_neighbors;}

    TRIANGLE_MESH& Triangle_Mesh() const
    {if(!triangle_mesh_up_to_date){Build_Mesh(triangle_mesh,0,cell_to_triangle_mapping,0);triangle_mesh_up_to_date=true;}
    return triangle_mesh;}

    ARRAY<int>& Cell_To_Triangle_Mapping() const
    {if(!triangle_mesh_up_to_date){Build_Mesh(triangle_mesh,0,cell_to_triangle_mapping,0);triangle_mesh_up_to_date=true;}
    return cell_to_triangle_mapping;}

    bool Outside(const VECTOR<T,2>& location) const
    {return uniform_grid.Outside(location);}

    template<class TV>
    TV Interpolate_Nodes(const ARRAY<TV>& u,const VECTOR<T,2>& X) const
    {Node_Locations();const RED_TRIANGLE<T>* red_leaf=Red_Leaf_Triangle(X);int i,j,k;
    if(red_leaf->Has_Green_Children()) red_leaf->green_children->elements(red_leaf->green_children->Green_Leaf_Triangle(X)).Get(i,j,k);
    else{i=red_leaf->Node(0);j=red_leaf->Node(1);k=red_leaf->Node(2);}
    VECTOR<T,3> weights=TRIANGLE_2D<T>::Barycentric_Coordinates(X,node_locations(i),node_locations(j),node_locations(k));
    return weights.x*u(i)+weights.y*u(j)+weights.z*u(k);}

    template<class TV>
    void Interpolate_Node_Values_To_New_Nodes(ARRAY<TV>& node_values,const int old_number_of_nodes)
    {for(int i=1;i<=uniform_grid.counts.x;i++)for(int j=1;j<=uniform_grid.counts.y;j++)if(elements(i,j))
        elements(i,j)->Interpolate_Node_Values_To_New_Nodes(node_values,old_number_of_nodes);}
    
    template<class TH>
    bool Refine_One_Level(TH* helper,bool (*refinement_criteria)(TH* helper,const RED_TRIANGLE<T>* triangle))
    {bool refinement=false;int old_number_of_cells=number_of_cells;
    for(int i=1;i<=uniform_grid.counts.x;i++)for(int j=1;j<=uniform_grid.counts.y;j++)if(elements(i,j))
        if(elements(i,j)->Refine(number_of_cells,number_of_nodes,old_number_of_cells,*this,helper,refinement_criteria))refinement=true;
    for(int i=1;i<=uniform_grid.counts.x;i++)for(int j=1;j<=uniform_grid.counts.y;j++)if(elements(i,j))
        elements(i,j)->Finalize_Green_Refinement(number_of_cells);
    return refinement;}

//#####################################################################
    void Initialize(const GRID<TV>& uniform_grid_input,const int maximum_depth_input);
    void Build_Mesh(TRIANGLE_MESH& triangle_mesh,const ARRAY<T>* phi,ARRAY<int>& cell_to_triangle_mapping,ARRAY<int>* node_to_particle_mapping) const;
    void Calculate_Node_Locations(ARRAY<VECTOR<T,2> >& node_locations) const;
    void Calculate_Node_Neighbors(ARRAY<ARRAY<int> >& node_neighbors) const;
    void Compact_Array_Indices(ARRAY<int>* cell_mapping_array,ARRAY<int>* node_mapping_array);
    const RED_TRIANGLE<T>* Red_Leaf_Triangle(const VECTOR<T,2>& location) const;
    const RED_TRIANGLE<T>* Root_Level_Triangle(const VECTOR<T,2>& location) const;
    void Create_Overlay_Grid(QUADTREE_GRID<T>& grid,const int number_of_ghost_cells,const bool use_nodes,const bool use_faces) const;
private:
    void Calculate_Cell_Pointer_From_Index(ARRAY<RED_TRIANGLE<T>*>& cell_pointer_from_index) const;
//#####################################################################
};
}
#endif
#endif
