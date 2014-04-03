//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license  contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BINTREE_CELL  
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __BINTREE_CELL__
#define __BINTREE_CELL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{ 

template<class T> class BINTREE_GRID;
template<class T> class BINTREE_CHILDREN;
template<class T1,class T2> class PAIR;

template<class T>
class BINTREE_CELL:public NONCOPYABLE
{
public:
    typedef VECTOR<int,1> TV_INT;
    enum CELL_NODES{LX=0,HX=1}; // 2 nodes of a cell (used for both children and nodes)
    enum AXIS{X_AXIS=0}; // used for recursive functions that depend on axis
    enum SIBLINGS{X_SIBLING=1}; // adding or subracting these values from an quadtree cell yields the neighbor in the given direction (if it's a sibling with the same parent)
    enum REFINE_ACTION{UNSPECIFIED=0,COARSEN=1,NONE=2,REFINE=3}; // none means that the cell should be left at the same level
    
    BINTREE_CHILDREN<T>* children; // pointer to the structure holding the 4 children of this cell - 0 if this is a leaf node
    BINTREE_CHILDREN<T>* owner; // pointer to the structure holding this cell and its 3 siblings
    int child_index; // the index of this child (gives its position in the BINTREE_CHILDREN structure, and its position relative to its siblings)
    int cell;
public:

    BINTREE_CELL()
    {}

    ~BINTREE_CELL()
    {if(Has_Children()) delete children;}

    void Initialize(BINTREE_CHILDREN<T>* owner_input,int& cell_input,const int child_index_input) 
    {assert(0<=child_index_input && child_index_input<2);owner=owner_input;children=0;child_index=child_index_input;cell=++cell_input;}
    
    static BINTREE_CELL<T>* Child(const int child,const BINTREE_CELL<T>* cell) 
    {return cell->Child(child);}

    static int Depth_Of_This_Cell(const BINTREE_CELL<T>* cell) 
    {return cell->Depth_Of_This_Cell();}

    bool Has_Children() const 
    {return children!=0;}

    int Nearest_Node(const VECTOR<T,1>& X) const
    {VECTOR<T,1> center=Center();return Node(X.x>center.x?1:0);}

    BINTREE_CELL<T>* Child(const int child) const
    {assert(children);return &children->children[child];}

    BINTREE_CELL<T>* Parent() const
    {return owner->parent;}

    int Depth_Of_This_Cell() const
    {return owner->childrens_depth;}

    int Maximum_Relative_Depth() const
    {if(!Has_Children()) return 1; int max_depth=0;for(int i=0;i<2;i++) max_depth=max(max_depth,Child(i)->Maximum_Relative_Depth()+1); return max_depth;}

    int Cell() const
    {return cell;}

    int& Cell()
    {return cell;}

    int& Node(const int node_index) const
    {return owner->Node(child_index,node_index);}

    int& Face(const int face_index) const
    {return owner->Face(child_index,face_index);}

    T Cell_Size() const
    {return DX().x;}

    T Face_Size() const
    {return 1;}

    const VECTOR<T,1>& DX() const
    {return owner->childrens_DX;}

    template<class TV> void Interpolate_Node_Values_To_Direct_Children(ARRAY<TV>& node_values)
    {assert(Has_Children()); node_values(children->Node(0,1))=(T).5*(node_values(Node(0))+node_values(Node(1)));}

    template<class TV> void Interpolate_Face_Values_To_Direct_Children(ARRAY<TV>& face_values)
    {assert(Has_Children()); face_values(children->Face(0,1))=(T).5*(face_values(Face(0))+face_values(Face(1)));}

//#####################################################################
    VECTOR<T,1> Center() const; // defined below for inlining
    const T Face_Diagonal_Length() const; // defined below for inlining
    const T Diagonal_Length() const; // defined below for inlining
    RANGE<VECTOR<T,1> > Bounding_Box() const; // defined below for inlining
//#####################################################################
    int Storage_Requirement() const;
    void Create_Cell_Compaction_Mapping(ARRAY<int>& mapping_array,int& cell_count);
    static BINTREE_CELL<T>* Create_Root(bool use_nodes,bool use_faces,int& number_of_cells,int& number_of_nodes,int& number_of_faces,const VECTOR<T,1>& root_center,const VECTOR<T,1>& root_DX);
    void Create_Children(int& number_of_cells,ARRAY<BINTREE_CELL<T>*>* new_cells,int& number_of_nodes,ARRAY<PAIR<BINTREE_CELL<T>*,int> >* new_nodes, 
        int& number_of_faces,ARRAY<PAIR<BINTREE_CELL<T>*,int> >* new_faces,const BINTREE_GRID<T>* grid);
    void Delete_Children();
    BINTREE_CELL<T>* Get_Neighbor(const int x_offset,const BINTREE_GRID<T>* grid,const ARRAY<VECTOR<BINTREE_CELL<T>*,2> >* neighbors=0,const bool same_depth_only=true) const;
    void Get_All_Face_Neighbors(int face_index,ARRAY<BINTREE_CELL<T>*>& face_neighbors,const BINTREE_GRID<T>* grid) const;
    template<class T_FACE_LOOKUP> void Interpolate_Face_Values_From_Direct_Children(ARRAY<typename T_FACE_LOOKUP::ELEMENT>& face_values,T_FACE_LOOKUP face_values_lookup);
private:
    static bool Can_Go_In_X_Direction(const int child_index,const int x_offset)
    {if(x_offset == 0) return true; else return (x_offset == 1 && child_index == 0) || (x_offset == -1 && child_index == 1);}
//#####################################################################
};
}
namespace PhysBAM{
//#####################################################################
// Function Center
//#####################################################################
template<class T> inline VECTOR<T,1> BINTREE_CELL<T>::
Center() const
{
    VECTOR<T,1>& DX=owner->childrens_DX;
    static const T center_lookup_table[]={-.5,.5};
    return VECTOR<T,1>(owner->parents_center.x+DX.x*center_lookup_table[child_index]);
}
//#####################################################################
// Function Face_Diagonal_Length
//#####################################################################
template<class T> inline const T BINTREE_CELL<T>::
Face_Diagonal_Length() const
{
    PHYSBAM_FATAL_ERROR("Face_Diagonal_Length() not implemented for BINTREE data structure");
    return owner->childrens_DX.x; // assumes square cells
}
//#####################################################################
// Function Diagonal_Length
//#####################################################################
template<class T> inline const T BINTREE_CELL<T>::
Diagonal_Length() const
{
    PHYSBAM_FATAL_ERROR("Diagonal_Length() not implemented for BINTREE data structure");
    return (T)root_two*owner->childrens_DX.x;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> inline RANGE<VECTOR<T,1> > BINTREE_CELL<T>::
Bounding_Box() const
{
    VECTOR<T,1> dx_over_two=DX()*(T).5,lower=Center()-dx_over_two,upper=Center()+dx_over_two;
    return RANGE<VECTOR<T,1> >(lower,upper);
}
//#####################################################################
}
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_CHILDREN.h>
#endif
#endif
