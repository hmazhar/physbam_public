//#####################################################################
// Copyright 2003-2004, Ron Fedkiw, Frank Losasso, Tamar Shinar, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license  contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OCTREE_CELL  
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __OCTREE_CELL__
#define __OCTREE_CELL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#ifdef HZ
#undef HZ
#endif
namespace PhysBAM{ 

template<class T> class OCTREE_GRID; 
template<class T> class OCTREE_CHILDREN;

template<class T>
class OCTREE_CELL:public NONCOPYABLE
{
    typedef VECTOR<int,3> TV_INT;
public:
    enum CELL_NODES{LX_LY_LZ=0,HX_LY_LZ=1,LX_HY_LZ=2,HX_HY_LZ=3,LX_LY_HZ=4,HX_LY_HZ=5,LX_HY_HZ=6,HX_HY_HZ=7}; // 8 nodes of a cell (used for both children and nodes)
    enum CELL_FACES{LX=0,HX=1,LY=2,HY=3,LZ=4,HZ=5}; // used to index the face structure
    enum AXIS{X_AXIS=0,Y_AXIS=1,Z_AXIS=2}; // used for recursive functions that depend on axis
    enum SIBLINGS{X_SIBLING=1,Y_SIBLING=2,Z_SIBLING=4}; // adding or subracting these values from an octree cell yields the neighbor in the given direction (if it's a sibling with the same parent)
    enum REFINE_ACTION{UNSPECIFIED=0,COARSEN=1,NONE=2,REFINE=3}; // none means that the cell should be left at the same level
    
    OCTREE_CHILDREN<T>* children; // pointer to the structure holding the 8 children of this cell - 0 if this is a leaf node
    OCTREE_CHILDREN<T>* owner; // pointer to the structure holding this cell and its 7 siblings
    unsigned char child_index; // the index of this child (gives its position in the OCTREE_CHILDREN structure, and its position relative to its siblings)
    int cell;
public:

    OCTREE_CELL()
    {}

    ~OCTREE_CELL();

    void Initialize(OCTREE_CHILDREN<T>* owner_input,int& number_of_cells,const unsigned char child_index_input) 
    {owner=owner_input;children=0;child_index=child_index_input;cell=++number_of_cells;}
    
    static OCTREE_CELL<T>* Child(const int child,const OCTREE_CELL<T>* cell) 
    {return cell->Child(child);}

    static int Depth_Of_This_Cell(const OCTREE_CELL<T>* cell) 
    {return cell->Depth_Of_This_Cell();}

    bool Has_Children() const 
    {return children!=0;}

    bool Intersection(const RANGE<VECTOR<T,3> >& box) const
    {VECTOR<T,3> center=Center();VECTOR<T,3> dx_over_two=(T).5*DX();
    RANGE<VECTOR<T,3> > cell_box(center.x-dx_over_two.x,center.x+dx_over_two.x,center.y-dx_over_two.y,center.y+dx_over_two.y,center.z-dx_over_two.z,center.z+dx_over_two.z);
    return cell_box.Lazy_Intersection(box);}

    int Nearest_Node(const VECTOR<T,3> &X) const
    {VECTOR<T,3> center=Center();return Node((X.x>center.x?1:0)|(X.y>center.y?2:0)|(X.z>center.z?4:0));}   

//#####################################################################
public:
    OCTREE_CELL<T>* Child(const int child) const; // defined below for inlining
    OCTREE_CELL<T>* Parent() const; // defined below for inlining
    int Depth_Of_This_Cell() const; // defined below for inlining
    int Maximum_Relative_Depth() const; // defined below for inlining
    int& Cell(); // defined below for inlining
    int Cell() const; // defined below for inlining
    int& Node(const int node_index) const; // defined below for inlining
    int& Face(const int face_index) const; // defined below for inlining
    VECTOR<T,3> Center() const; // defined below for inlining
    const VECTOR<T,3>& DX() const; // defined below for inlining
    const T Face_Diagonal_Length() const; // defined below for inlining
    const T Diagonal_Length() const; // defined below for inlining
    RANGE<VECTOR<T,3> > Bounding_Box() const; // defined below for inlining
    T Cell_Size() const; // defined below for inlining
    T Face_Size() const; // defined below for inlining
public:
    void Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count);
public:
    int Storage_Requirement() const;
    static OCTREE_CELL<T>* Create_Root(bool use_nodes,bool use_faces,int& number_of_cells,int& number_of_nodes,int& number_of_faces,const VECTOR<T,3>& root_center,const VECTOR<T,3>& root_DX);
    void Create_Children(int& number_of_cells,ARRAY<OCTREE_CELL<T>*>* new_cells,int& number_of_nodes,ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_nodes, 
        int& number_of_faces,ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_faces,const OCTREE_GRID<T>* grid);
    void Delete_Children();
    OCTREE_CELL<T>* Get_Neighbor(const int x_offset,const int y_offset,const int z_offset,const OCTREE_GRID<T>* grid,const ARRAY<VECTOR<OCTREE_CELL<T>*,6> >* neighbors=0,const bool same_depth_only=true) const;
    void Get_All_Face_Neighbors(const int face_index,ARRAY<OCTREE_CELL<T>*>& face_neighbors,const OCTREE_GRID<T>* grid) const;
    template<class TV> void Interpolate_Node_Values_To_Direct_Children(ARRAY<TV>& node_values);
    template<class TV> void Interpolate_Face_Values_To_Direct_Children(ARRAY<TV>& face_values);
private:
    static bool Can_Go_In_X_Direction(const int child_index,const int x_offset);
    static bool Can_Go_In_Y_Direction(const int child_index,const int y_offset);
    static bool Can_Go_In_Z_Direction(const int child_index,const int z_offset);
//#####################################################################
};
}
namespace PhysBAM{
//#####################################################################
// Function ~OCTREE_CELL
//#####################################################################
template<class T> inline OCTREE_CELL<T>::
~OCTREE_CELL()
{
    if(Has_Children()) delete children;
}
//#####################################################################
// Function Child
//#####################################################################
template<class T> inline OCTREE_CELL<T>* OCTREE_CELL<T>::
Child(const int child) const 
{
    assert(children);
    return &children->children[child];
}
//#####################################################################
// Function Parent
//#####################################################################
template<class T> inline OCTREE_CELL<T>* OCTREE_CELL<T>::
Parent() const 
{
    return owner->parent;
}
//#####################################################################
// Function Depth_Of_This_Cell
//#####################################################################
template<class T> inline int OCTREE_CELL<T>::
Depth_Of_This_Cell() const
{
    return owner->childrens_depth;
}
//#####################################################################
// Function Maximum_Relative_Depth
//#####################################################################
template<class T> inline int OCTREE_CELL<T>::
Maximum_Relative_Depth() const
{
    if(!Has_Children()) return 1;
    int max_depth=0;for(int i=0;i<8;i++) max_depth=max(max_depth,Child(i)->Maximum_Relative_Depth()+1);
    return max_depth;
}
//#####################################################################
// Function Center
//#####################################################################
template<class T> inline VECTOR<T,3> OCTREE_CELL<T>::
Center() const
{
    VECTOR<T,3>& DX=owner->childrens_DX;
    static const T center_lookup_table[][3]={{-.5,-.5,-.5},{.5,-.5,-.5},{-.5,.5,-.5},{.5,.5,-.5},{-.5,-.5,.5},{.5,-.5,.5},{-.5,.5,.5},{.5,.5,.5}};
    return VECTOR<T,3>(owner->parents_center.x+DX.x*center_lookup_table[child_index][0],owner->parents_center.y+DX.y*center_lookup_table[child_index][1],owner->parents_center.z+DX.z*center_lookup_table[child_index][2]);
}
//#####################################################################
// Function DX
//#####################################################################
template<class T> inline const VECTOR<T,3>& OCTREE_CELL<T>::
DX() const
{
    return owner->childrens_DX;
}
//#####################################################################
// Function Diagonal_Length
//#####################################################################
template<class T> inline const T OCTREE_CELL<T>::
Face_Diagonal_Length() const
{
    return (T)root_two*owner->childrens_DX.x; // assumes square cells
}
//#####################################################################
// Function Diagonal_Length
//#####################################################################
template<class T> inline const T OCTREE_CELL<T>::
Diagonal_Length() const
{
    return (T)root_three*owner->childrens_DX.x; // assumes square cells
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> inline RANGE<VECTOR<T,3> > OCTREE_CELL<T>::
Bounding_Box() const
{
    VECTOR<T,3>& DX=owner->childrens_DX;
    static const T minimum_corner_lookup_table[][3]={{-1,-1,-1},{0,-1,-1},{-1,0,-1},{0,0,-1},{-1,-1,0},{0,-1,0},{-1,0,0},{0,0,0}};
    VECTOR<T,3> minimum_corner(owner->parents_center.x+DX.x*minimum_corner_lookup_table[child_index][0],owner->parents_center.y+DX.y*minimum_corner_lookup_table[child_index][1],owner->parents_center.z+DX.z*minimum_corner_lookup_table[child_index][2]);
    return RANGE<VECTOR<T,3> >(minimum_corner,minimum_corner+DX);
}
//#####################################################################
// Function Cell_Size
//#####################################################################
template<class T> inline T OCTREE_CELL<T>::
Cell_Size() const
{
    return DX().x*DX().y*DX().z;
}
//#####################################################################
// Function Face_Size
//#####################################################################
template<class T> inline T OCTREE_CELL<T>::
Face_Size() const
{
    return DX().x*DX().x;
}
//#####################################################################
// Function Cell
//#####################################################################
template<class T> inline int& OCTREE_CELL<T>::
Cell()
{
    return cell;
}
//#####################################################################
// Function Cell
//#####################################################################
template<class T> inline int OCTREE_CELL<T>::
Cell() const
{
    return cell;
}
//#####################################################################
// Function Node
//#####################################################################
template<class T> inline int& OCTREE_CELL<T>::
Node(const int node_index) const 
{
    return owner->Node(child_index,node_index);
}
//#####################################################################
// Function Face
//#####################################################################
template<class T> inline int& OCTREE_CELL<T>::
Face(const int face_index) const 
{
    return owner->Face(child_index,face_index);
}
//#####################################################################
}
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CHILDREN.h>
#endif
#endif
