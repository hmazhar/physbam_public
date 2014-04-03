//#####################################################################
// Copyright 2004-2005, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license  contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUADTREE_CELL  
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __QUADTREE_CELL__
#define __QUADTREE_CELL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{ 

template<class T> class QUADTREE_GRID;
template<class T> class QUADTREE_CHILDREN;
template<class T1,class T2> class PAIR;

template<class T>
class QUADTREE_CELL:public NONCOPYABLE
{
public:
    typedef VECTOR<int,2> TV_INT;
    enum CELL_NODES{LX_LY=0,HX_LY=1,LX_HY=2,HX_HY=3}; // 4 nodes of a cell (used for both children and nodes)
    enum CELL_FACES{LX=0,HX=1,LY=2,HY=3}; // used to index the face structure
    enum AXIS{X_AXIS=0,Y_AXIS=1}; // used for recursive functions that depend on axis
    enum SIBLINGS{X_SIBLING=1,Y_SIBLING=2}; // adding or subracting these values from an quadtree cell yields the neighbor in the given direction (if it's a sibling with the same parent)
    enum REFINE_ACTION{UNSPECIFIED=0,COARSEN=1,NONE=2,REFINE=3}; // none means that the cell should be left at the same level
    
    QUADTREE_CHILDREN<T>* children; // pointer to the structure holding the 4 children of this cell - 0 if this is a leaf node
    QUADTREE_CHILDREN<T>* owner; // pointer to the structure holding this cell and its 3 siblings
    unsigned char child_index; // the index of this child (gives its position in the QUADTREE_CHILDREN structure, and its position relative to its siblings)
    int cell;
public:

    QUADTREE_CELL()
    {}

    ~QUADTREE_CELL();

    void Initialize(QUADTREE_CHILDREN<T>* owner_input,int& cell_input,const unsigned char child_index_input) 
    {owner=owner_input;children=0;child_index=child_index_input;cell=++cell_input;}
    
    static QUADTREE_CELL<T>* Child(const int child,const QUADTREE_CELL<T>* cell) 
    {return cell->Child(child);}

    static int Depth_Of_This_Cell(const QUADTREE_CELL<T>* cell) 
    {return cell->Depth_Of_This_Cell();}

    bool Has_Children() const 
    {return children!=0;}

    int Nearest_Node(const VECTOR<T,2>& X) const
    {VECTOR<T,2> center=Center();return Node((X.x>center.x?1:0)|(X.y>center.y?2:0));}

//#####################################################################
public:
    QUADTREE_CELL<T>* Child(const int child) const; // defined below for inlining
    QUADTREE_CELL<T>* Parent() const; // defined below for inlining
    int Depth_Of_This_Cell() const; // defined below for inlining
    int Maximum_Relative_Depth() const; // defined below for inlining
    int Cell() const; // defined below for inlining
    int& Cell(); // defined below for inlining
    int& Node(const int node_index) const; // defined below for inlining
    int& Face(const int face_index) const; // defined below for inlining
    VECTOR<T,2> Center() const; // defined below for inlining
    const VECTOR<T,2>& DX() const; // defined below for inlining
    const T Face_Diagonal_Length() const; // defined below for inlining
    const T Diagonal_Length() const; // defined below for inlining
    RANGE<VECTOR<T,2> > Bounding_Box() const; // defined below for inlining
    T Cell_Size() const; // defined below for inlining
    T Face_Size() const; // defined below for inlining
public:
    void Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count);
public:
    int Storage_Requirement() const;
    static QUADTREE_CELL<T>* Create_Root(bool use_nodes,bool use_faces,int& number_of_cells,int& number_of_nodes,int& number_of_faces,const VECTOR<T,2>& root_center,const VECTOR<T,2>& root_DX);
    void Create_Children(int& number_of_cells,ARRAY<QUADTREE_CELL<T>*>* new_cells,int& number_of_nodes,ARRAY<PAIR<QUADTREE_CELL<T>*,int> >* new_nodes, 
        int& number_of_faces,ARRAY<PAIR<QUADTREE_CELL<T>*,int> >* new_faces,const QUADTREE_GRID<T>* grid);
    void Delete_Children();
    QUADTREE_CELL<T>* Get_Neighbor(const int x_offset,const int y_offset,const QUADTREE_GRID<T>* grid,const ARRAY<VECTOR<QUADTREE_CELL<T>*,4> >* neighbors=0,const bool same_depth_only=true) const;
    void Get_All_Face_Neighbors(int face_index,ARRAY<QUADTREE_CELL<T>*>& face_neighbors,const QUADTREE_GRID<T>* grid) const;
    template<class TV> void Interpolate_Node_Values_To_Direct_Children(ARRAY<TV>& node_values);
    template<class TV> void Interpolate_Face_Values_To_Direct_Children(ARRAY<TV>& face_values);
    template<class T_FACE_LOOKUP> void Interpolate_Face_Values_From_Direct_Children(ARRAY<typename T_FACE_LOOKUP::ELEMENT>& face_values,T_FACE_LOOKUP face_values_lookup);
private:
    static bool Can_Go_In_X_Direction(const int child_index,const int x_offset);
    static bool Can_Go_In_Y_Direction(const int child_index,const int y_offset);
//#####################################################################
};
}
namespace PhysBAM{
//#####################################################################
// Function ~QUADTREE_CELL
//#####################################################################
template<class T> inline QUADTREE_CELL<T>::
~QUADTREE_CELL()
{
    if(Has_Children()) delete children;
}
//#####################################################################
// Function Child
//#####################################################################
template<class T> inline QUADTREE_CELL<T>* QUADTREE_CELL<T>::
Child(const int child) const 
{
    assert(children);
    return &children->children[child];
}
//#####################################################################
// Function Parent
//#####################################################################
template<class T> inline QUADTREE_CELL<T>* QUADTREE_CELL<T>::
Parent() const 
{
    return owner->parent;
}
//#####################################################################
// Function Depth_Of_This_Cell
//#####################################################################
template<class T> inline int QUADTREE_CELL<T>::
Depth_Of_This_Cell() const
{
    return owner->childrens_depth;
}
//#####################################################################
// Function Maximum_Relative_Depth
//#####################################################################
template<class T> inline int QUADTREE_CELL<T>::
Maximum_Relative_Depth() const
{
    if(!Has_Children()) return 1;
    int max_depth=0;for(int i=0;i<4;i++) max_depth=max(max_depth,Child(i)->Maximum_Relative_Depth()+1);
    return max_depth;
}
//#####################################################################
// Function Center
//#####################################################################
template<class T> inline VECTOR<T,2> QUADTREE_CELL<T>::
Center() const
{
    VECTOR<T,2>& DX=owner->childrens_DX;
    static const T center_lookup_table[][3]={{-.5,-.5},{.5,-.5},{-.5,.5},{.5,.5}};
    return VECTOR<T,2>(owner->parents_center.x+DX.x*center_lookup_table[child_index][0],owner->parents_center.y+DX.y*center_lookup_table[child_index][1]);
}
//#####################################################################
// Function DX
//#####################################################################
template<class T> inline const VECTOR<T,2>& QUADTREE_CELL<T>::
DX() const
{
    return owner->childrens_DX;
}
//#####################################################################
// Function Face_Diagonal_Length
//#####################################################################
template<class T> inline const T QUADTREE_CELL<T>::
Face_Diagonal_Length() const
{
    return owner->childrens_DX.x; // assumes square cells
}
//#####################################################################
// Function Diagonal_Length
//#####################################################################
template<class T> inline const T QUADTREE_CELL<T>::
Diagonal_Length() const
{
    return (T)root_two*owner->childrens_DX.x;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> inline RANGE<VECTOR<T,2> > QUADTREE_CELL<T>::
Bounding_Box() const
{
    VECTOR<T,2> dx_over_two=DX()*(T).5,lower=Center()-dx_over_two,upper=Center()+dx_over_two;
    return RANGE<VECTOR<T,2> >(lower,upper);
}
//#####################################################################
// Function Cell_Size
//#####################################################################
template<class T> inline T QUADTREE_CELL<T>::
Cell_Size() const
{
    return DX().x*DX().y;
}
//#####################################################################
// Function Face_Size
//#####################################################################
template<class T> inline T QUADTREE_CELL<T>::
Face_Size() const
{
    return DX().x;
}
//#####################################################################
// Function Cell
//#####################################################################
template<class T> inline int& QUADTREE_CELL<T>::
Cell()
{
    return cell;
}
//#####################################################################
// Function Cell
//#####################################################################
template<class T> inline int QUADTREE_CELL<T>::
Cell() const
{
    return cell;
}
//#####################################################################
// Function Node
//#####################################################################
template<class T> inline int& QUADTREE_CELL<T>::
Node(const int node_index) const 
{
    return owner->Node(child_index,node_index);
}
//#####################################################################
// Function Face
//#####################################################################
template<class T> inline int& QUADTREE_CELL<T>::
Face(const int face_index) const 
{
    return owner->Face(child_index,face_index);
}
//#####################################################################
}
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CHILDREN.h>
#endif
#endif
