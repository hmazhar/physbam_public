//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license  contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_TETRAHEDRON  
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __RED_TETRAHEDRON__
#define __RED_TETRAHEDRON__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Red_Green/TETRAHEDRAL_GROUP.h>
namespace PhysBAM{ 

class TETRAHEDRON_MESH;
template<class T> class RED_CHILDREN_3D;
template<class T> class GREEN_CHILDREN_3D;
template<class T> class RED_GREEN_GRID_3D;

template<class T>
class RED_TETRAHEDRON:public NONCOPYABLE
{
public:
    RED_CHILDREN_3D<T>* red_children; // pointer to the structure holding the 8 red children of this tetrahedron or 0 if no red children
    GREEN_CHILDREN_3D<T>* green_children; // pointer to the structure holding the green children of this tetrahedron or 0 if leaf or red refined
    RED_CHILDREN_3D<T>* owner; // pointer to the structure holding this tetrahedron and its 7 siblings
    unsigned char child_index; // the index of this child (gives its position in the RED_CHILDREN_3D structure, and its position relative to its siblings)
    int cell;

    RED_TETRAHEDRON()
    {}

    ~RED_TETRAHEDRON();
    
    void Initialize(RED_CHILDREN_3D<T>* owner_input)
    {owner=owner_input;red_children=0;green_children=0;child_index=0;cell=0;}
    
    void Initialize(RED_CHILDREN_3D<T>* owner_input,int& cell_input,const unsigned char child_index_input) 
    {owner=owner_input;red_children=0;green_children=0;child_index=child_index_input;cell=++cell_input;}
        
    bool Has_Red_Children() const
    {return red_children!=0;}

    bool Has_Green_Children() const
    {return green_children!=0;}

    int Maximum_Relative_Depth() const
    {if(!Has_Red_Children()) return 1+Has_Green_Children();
    int max_depth=0;for(int i=0;i<8;i++) max_depth=max(max_depth,Red_Child(i)->Maximum_Relative_Depth()+1);
    return max_depth;}

    template<class TV>
    void Interpolate_Node_Values_To_Direct_Red_Children(ARRAY<TV>& node_values,const int old_number_of_nodes=0)
    {assert(Has_Red_Children());
    if(red_children->Midpoint(0)>old_number_of_nodes)node_values(red_children->Midpoint(0))=(T).5*(node_values(Node(0))+node_values(Node(1)));
    if(red_children->Midpoint(1)>old_number_of_nodes)node_values(red_children->Midpoint(1))=(T).5*(node_values(Node(2))+node_values(Node(3)));
    if(red_children->Midpoint(2)>old_number_of_nodes)node_values(red_children->Midpoint(2))=(T).5*(node_values(Node(0))+node_values(Node(2)));
    if(red_children->Midpoint(3)>old_number_of_nodes)node_values(red_children->Midpoint(3))=(T).5*(node_values(Node(0))+node_values(Node(3)));
    if(red_children->Midpoint(4)>old_number_of_nodes)node_values(red_children->Midpoint(4))=(T).5*(node_values(Node(1))+node_values(Node(3)));
    if(red_children->Midpoint(5)>old_number_of_nodes)node_values(red_children->Midpoint(5))=(T).5*(node_values(Node(1))+node_values(Node(2)));}

    template<class TV>
    void Interpolate_Node_Values_To_Direct_Green_Children(ARRAY<TV>& node_values,const int old_number_of_nodes=0)
    {assert(Has_Green_Children());
    if(green_children->midpoints[0]>old_number_of_nodes)node_values(green_children->midpoints[0])=(T).5*(node_values(Node(0))+node_values(Node(1)));
    if(green_children->midpoints[1]>old_number_of_nodes)node_values(green_children->midpoints[1])=(T).5*(node_values(Node(2))+node_values(Node(3)));
    if(green_children->midpoints[2]>old_number_of_nodes)node_values(green_children->midpoints[2])=(T).5*(node_values(Node(0))+node_values(Node(2)));
    if(green_children->midpoints[3]>old_number_of_nodes)node_values(green_children->midpoints[3])=(T).5*(node_values(Node(0))+node_values(Node(3)));
    if(green_children->midpoints[4]>old_number_of_nodes)node_values(green_children->midpoints[4])=(T).5*(node_values(Node(1))+node_values(Node(3)));
    if(green_children->midpoints[5]>old_number_of_nodes)node_values(green_children->midpoints[5])=(T).5*(node_values(Node(1))+node_values(Node(2)));}

    template<class TV>
    void Interpolate_Node_Values_To_All_Children(ARRAY<TV>& node_values)
    {if(Has_Red_Children()){Interpolate_Node_Values_To_Direct_Red_Children(node_values);
        for(int i=0;i<8;i++)Red_Child(i)->Interpolate_Node_Values_To_All_Children(node_values);}
    else if(Has_Green_Children())Interpolate_Node_Values_To_Direct_Green_Children(node_values);}
    
    template<class TV>
    void Interpolate_Node_Values_To_New_Nodes(ARRAY<TV>& node_values,const int old_number_of_nodes)
    {if(Has_Red_Children()){Interpolate_Node_Values_To_Direct_Red_Children(node_values,old_number_of_nodes);
        for(int i=0;i<8;i++)Red_Child(i)->Interpolate_Node_Values_To_New_Nodes(node_values,old_number_of_nodes);}
    else if(Has_Green_Children())Interpolate_Node_Values_To_Direct_Green_Children(node_values,old_number_of_nodes);}
    
    template<class TH>
    bool Refine(int& number_of_cells,int& number_of_nodes,const int old_number_of_cells,const RED_GREEN_GRID_3D<T>& grid,TH* helper,
        bool (*refinement_criteria)(TH* helper,const RED_TETRAHEDRON<T>* tetrahedron))
    {if(cell>old_number_of_cells)return false;// make sure that refinement happens one level at a time
    if(Has_Red_Children()){bool refinement=false;
        for(int i=0;i<8;i++)if(Red_Child(i)->Refine(number_of_cells,number_of_nodes,old_number_of_cells,grid,helper,refinement_criteria))refinement=true;
        return refinement;}
    if(!refinement_criteria(helper,this))return false;
    Red_Refine(number_of_cells,number_of_nodes,grid);return true;}

    const RED_TETRAHEDRON<T>* Red_Leaf_Tetrahedron(const VECTOR<T,3>& location) const
    {return Red_Descendant(location,1000);} // go down as far as possible

//#####################################################################
    void Clean_Memory(); // defined below for inlining
    RED_TETRAHEDRON<T>* Red_Child(const int child) const; // defined below for inlining
    RED_TETRAHEDRON<T>* Parent() const; // defined below for inlining
    int Depth() const; // defined below for inlining
    int& Node(const int node_index) const; // defined below for inlining
    VECTOR<T,3> Center() const; // defined below for inlining
    TETRAHEDRAL_GROUP<T> Orientation() const; // defined below for inlining
    T Size() const; // defined below for inlining
    void Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count);
    void Red_Refine(int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_3D<T>& grid);
    void Propogate_Refinement(const int midpoint,int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_3D<T>& grid);
    void Simulate_Green_Refine(int& number_of_cells,int& number_of_nodes,const RED_GREEN_GRID_3D<T>& grid);
    void Finalize_Green_Refinement(int& number_of_cells);
    void Delete_Children();
    RED_TETRAHEDRON<T>* Get_Red_Neighbor(const int face,const RED_GREEN_GRID_3D<T>& grid,const bool same_depth_only=true) const;
    RED_TETRAHEDRON<T>* Get_Red_Neighbor_Ancestor(const int face,const RED_GREEN_GRID_3D<T>& grid) const;
    void Get_Red_Edge_Neighbors(const int edge,ARRAY<RED_TETRAHEDRON<T>*>& neighbors,const RED_GREEN_GRID_3D<T>& grid,const bool same_depth_only);
    int Get_Neighbor_Midpoint(const int midpoint,const RED_GREEN_GRID_3D<T>& grid);
    void Build_Tetrahedron_Mesh(TETRAHEDRON_MESH& tetrahedron_mesh,const ARRAY<T>* phi,ARRAY<int>& cell_to_tetrahedron_mapping,ARRAY<int>* node_to_particle_mapping) const;
    VECTOR<T,3> Node_Location(const int node_index) const;
    VECTOR<T,3> Midpoint_Location(const int midpoint) const;
    bool Outside(const VECTOR<T,3>& location,const T thickness_over_2=0) const;
    const RED_TETRAHEDRON<T>* Red_Descendant(const VECTOR<T,3>& location,const int depth) const;
    RED_TETRAHEDRON<T>* Red_Descendant(const VECTOR<T,3>& location,const int depth);
private:
    int Closest_Child(const VECTOR<T,3>& location) const;
    int Midpoint_From_Location(const VECTOR<T,3>& location) const;
    int Face_From_Face_Jump(const VECTOR<int,3>& face_jump) const;
//#####################################################################
};
}
namespace PhysBAM{
//#####################################################################
// Function ~RED_TETRAHEDRON
//#####################################################################
template<class T> inline RED_TETRAHEDRON<T>::
~RED_TETRAHEDRON()
{
    delete red_children;delete green_children;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T> inline void RED_TETRAHEDRON<T>::
Clean_Memory()
{
    if(red_children)red_children->Clean_Memory();delete red_children;red_children=0;delete green_children;green_children=0;
}
//#####################################################################
// Function Child
//#####################################################################
template<class T> inline RED_TETRAHEDRON<T>* RED_TETRAHEDRON<T>::
Red_Child(const int child) const 
{
    assert(red_children);
    return &red_children->children[child];
}
//#####################################################################
// Function Parent
//#####################################################################
template<class T> inline RED_TETRAHEDRON<T>* RED_TETRAHEDRON<T>::
Parent() const 
{
    return owner->parent;
}
//#####################################################################
// Function Depth_Of_This_Cell
//#####################################################################
template<class T> inline int RED_TETRAHEDRON<T>::
Depth() const
{
    return owner->childrens_depth;
}
//#####################################################################
// Function Node
//#####################################################################
template<class T> inline int& RED_TETRAHEDRON<T>::
Node(const int node_index) const 
{
    return owner->Node(child_index,node_index);
}
//#####################################################################
// Function Center
//#####################################################################
template<class T> inline VECTOR<T,3> RED_TETRAHEDRON<T>::
Center() const 
{
    static const VECTOR<T,3> center_lookup_table[8]={VECTOR<T,3>(-2,-1,0),VECTOR<T,3>(2,-1,0),VECTOR<T,3>(0,1,-2),VECTOR<T,3>(0,1,2),
        VECTOR<T,3>(-1,0,0),VECTOR<T,3>(0,0,1),VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,-1)};
    return owner->center+owner->childrens_size*owner->orientation.Rotate(center_lookup_table[child_index]);
}
//#####################################################################
// Function Orientation
//#####################################################################
template<class T> inline TETRAHEDRAL_GROUP<T> RED_TETRAHEDRON<T>::
Orientation() const 
{
    return owner->Child_Orientation(child_index);
}
//#####################################################################
// Function Size
//#####################################################################
template<class T> inline T RED_TETRAHEDRON<T>::
Size() const 
{
    return owner->childrens_size;
}
//#####################################################################
}
#include <PhysBAM_Geometry/Red_Green/GREEN_CHILDREN_3D.h>
#include <PhysBAM_Geometry/Red_Green/RED_CHILDREN_3D.h>
#endif
#endif
