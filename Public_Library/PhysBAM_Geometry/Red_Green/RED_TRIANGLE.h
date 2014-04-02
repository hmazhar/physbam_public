//#####################################################################
// Copyright 2004, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license  contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_TRIANGLE  
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __RED_TRIANGLE__
#define __RED_TRIANGLE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

class TRIANGLE_MESH;
template<class T> class RED_CHILDREN_2D;
template<class T> class GREEN_CHILDREN_2D;
template<class T> class RED_GREEN_GRID_2D;
template<class T> class QUADTREE_GRID;
template<class T> class QUADTREE_CELL;

template<class T>
class RED_TRIANGLE:public NONCOPYABLE
{
public:
    RED_CHILDREN_2D<T>* red_children; // pointer to the structure holding the 4 red children of this triangle or 0 if no red children
    GREEN_CHILDREN_2D<T>* green_children; // pointer to the structure holding the green children of this triangle or 0 if leaf or red refined
    RED_CHILDREN_2D<T>* owner; // pointer to the structure holding this triangle and its 3 siblings
    unsigned char child_index; // the index of this child (gives its position in the RED_CHILDREN_2D structure, and its position relative to its siblings)
    int cell;

    RED_TRIANGLE()
    {}

    ~RED_TRIANGLE();
    
    void Initialize(RED_CHILDREN_2D<T>* owner_input) 
    {owner=owner_input;red_children=0;green_children=0;child_index=0;cell=0;}
    
    void Initialize(RED_CHILDREN_2D<T>* owner_input,int& cell_input,const unsigned char child_index_input) 
    {owner=owner_input;red_children=0;green_children=0;child_index=child_index_input;cell=++cell_input;}
        
    bool Has_Red_Children() const
    {return red_children!=0;}

    bool Has_Green_Children() const
    {return green_children!=0;}

    int Maximum_Relative_Depth() const
    {if(!Has_Red_Children()) return 1+Has_Green_Children();
    int max_depth=0;for(int i=0;i<4;i++) max_depth=max(max_depth,Red_Child(i)->Maximum_Relative_Depth()+1);
    return max_depth;}

    template<class TV>
    void Interpolate_Node_Values_To_Direct_Red_Children(ARRAY<TV>& node_values,const int old_number_of_nodes=0)
    {assert(Has_Red_Children());
    if(red_children->Midpoint(0)>old_number_of_nodes)node_values(red_children->Midpoint(0))=(T).5*(node_values(Node(0))+node_values(Node(1)));
    if(red_children->Midpoint(1)>old_number_of_nodes)node_values(red_children->Midpoint(1))=(T).5*(node_values(Node(1))+node_values(Node(2)));
    if(red_children->Midpoint(2)>old_number_of_nodes)node_values(red_children->Midpoint(2))=(T).5*(node_values(Node(2))+node_values(Node(0)));}

    template<class TV>
    void Interpolate_Node_Values_To_Direct_Green_Children(ARRAY<TV>& node_values,const int old_number_of_nodes=0)
    {assert(Has_Green_Children());
    if(green_children->midpoints[0]>old_number_of_nodes)node_values(green_children->midpoints[0])=(T).5*(node_values(Node(0))+node_values(Node(1)));
    if(green_children->midpoints[1]>old_number_of_nodes)node_values(green_children->midpoints[1])=(T).5*(node_values(Node(1))+node_values(Node(2)));
    if(green_children->midpoints[2]>old_number_of_nodes)node_values(green_children->midpoints[2])=(T).5*(node_values(Node(2))+node_values(Node(0)));}

    template<class TV>
    void Interpolate_Node_Values_To_All_Children(ARRAY<TV>& node_values)
    {if(Has_Red_Children()){Interpolate_Node_Values_To_Direct_Red_Children(node_values);
        for(int i=0;i<4;i++)Red_Child(i)->Interpolate_Node_Values_To_All_Children(node_values);}
    else if(Has_Green_Children())Interpolate_Node_Values_To_Direct_Green_Children(node_values);}
    
    template<class TV>
    void Interpolate_Node_Values_To_New_Nodes(ARRAY<TV>& node_values,const int old_number_of_nodes)
    {if(Has_Red_Children()){Interpolate_Node_Values_To_Direct_Red_Children(node_values,old_number_of_nodes);
        for(int i=0;i<4;i++)Red_Child(i)->Interpolate_Node_Values_To_New_Nodes(node_values,old_number_of_nodes);}
    else if(Has_Green_Children())Interpolate_Node_Values_To_Direct_Green_Children(node_values,old_number_of_nodes);}
    
    template<class TH>
    bool Refine(int& number_of_cells,int& number_of_nodes,const int old_number_of_cells,const RED_GREEN_GRID_2D<T>& grid,TH* helper,
        bool (*refinement_criteria)(TH* helper,const RED_TRIANGLE<T>* triangle))
    {if(cell>old_number_of_cells)return false;// make sure that refinement happens one level at a time
    if(Has_Red_Children()){bool refinement=false;
        for(int i=0;i<4;i++)if(Red_Child(i)->Refine(number_of_cells,number_of_nodes,old_number_of_cells,grid,helper,refinement_criteria))refinement=true;
        return refinement;}
    if(!refinement_criteria(helper,this))return false;
    Red_Refine(number_of_cells,number_of_nodes,grid);return true;}

//#####################################################################
    void Clean_Memory(); // defined below for inlining
    RED_TRIANGLE<T>* Red_Child(const int child) const; // defined below for inlining
    RED_TRIANGLE<T>* Parent() const; // defined below for inlining
    int Depth() const; // defined below for inlining
    int& Node(const int node_index) const; // defined below for inlining
    VECTOR<T,2> Center() const; // defined below for inlining
    int Orientation() const; // defined below for inlining
    void Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count);
    void Red_Refine(int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_2D<T>& grid);
    void Propogate_Refinement(const int midpoint,int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_2D<T>& grid);
    void Simulate_Green_Refine(int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_2D<T>& grid);
    void Finalize_Green_Refinement(int& number_of_cells);
    void Delete_Children();
    RED_TRIANGLE<T>* Get_Red_Neighbor(const int edge,const RED_GREEN_GRID_2D<T>& grid,const bool same_depth_only=true) const;
    int Get_Neighbor_Midpoint(const int midpoint,const RED_GREEN_GRID_2D<T>& grid) const;
    void Build_Triangle_Mesh(TRIANGLE_MESH& triangle_mesh,const ARRAY<T>* phi,ARRAY<int>& cell_to_triangle_mapping,ARRAY<int>* node_to_particle_mapping) const;
    VECTOR<T,2> Node_Location(int node_index) const;
    bool Outside(const VECTOR<T,2>& location,const T thickness_over_2=0) const;
    const RED_TRIANGLE<T>* Red_Leaf_Triangle(const VECTOR<T,2>& location) const;
    void Calculate_Cell_Pointer_From_Index(ARRAY<RED_TRIANGLE<T>*>& cell_pointer_from_index);
    void Add_Ordered_Neighbors(ARRAY<ARRAY<int> >& neighbors,ARRAY<ARRAY<int> >& neighbor_links) const;
    void Create_Overlay_Grid(QUADTREE_GRID<T>& overlay_grid,QUADTREE_CELL<T>* cell1,QUADTREE_CELL<T>* cell2) const;
//#####################################################################
};
}
namespace PhysBAM{
//#####################################################################
// Function ~RED_TRIANGLE
//#####################################################################
template<class T> inline RED_TRIANGLE<T>::
~RED_TRIANGLE()
{
    delete red_children;delete green_children;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T> inline void RED_TRIANGLE<T>::
Clean_Memory()
{
    if(red_children)red_children->Clean_Memory();delete red_children;red_children=0;delete green_children;green_children=0;
}
//#####################################################################
// Function Child
//#####################################################################
template<class T> inline RED_TRIANGLE<T>* RED_TRIANGLE<T>::
Red_Child(const int child) const 
{
    assert(red_children);
    return &red_children->children[child];
}
//#####################################################################
// Function Parent
//#####################################################################
template<class T> inline RED_TRIANGLE<T>* RED_TRIANGLE<T>::
Parent() const 
{
    return owner->parent;
}
//#####################################################################
// Function Depth_Of_This_Cell
//#####################################################################
template<class T> inline int RED_TRIANGLE<T>::
Depth() const
{
    return owner->childrens_depth;
}
//#####################################################################
// Function Node
//#####################################################################
template<class T> inline int& RED_TRIANGLE<T>::
Node(const int node_index) const 
{
    return owner->Node(child_index,node_index);
}
//#####################################################################
// Function Center
//#####################################################################
template<class T> inline VECTOR<T,2> RED_TRIANGLE<T>::
Center() const 
{
    static const VECTOR<T,2> center_lookup_table[4]={VECTOR<T,2>(-2,-1),VECTOR<T,2>(2,-1),VECTOR<T,2>(1,0),VECTOR<T,2>(-1,0)};
    return owner->center+owner->childrens_size*center_lookup_table[child_index].Rotate_Counterclockwise_Multiple_90(owner->orientation);
}
//#####################################################################
// Function Orientation
//#####################################################################
template<class T> inline int RED_TRIANGLE<T>::
Orientation() const 
{
    return owner->Child_Orientation(child_index);
}
//#####################################################################
}
#include <PhysBAM_Geometry/Red_Green/GREEN_CHILDREN_2D.h>
#include <PhysBAM_Geometry/Red_Green/RED_CHILDREN_2D.h>
#endif
#endif
