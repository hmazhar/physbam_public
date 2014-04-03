//#####################################################################
// Copyright 2004-2009, Zhaosheng Bao, Ron Fedkiw, Jon Gretarsson, Geoffrey Irving, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOX_HIERARCHY
//##################################################################### 
#ifndef __BOX_HIERARCHY__
#define __BOX_HIERARCHY__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/ZERO.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
namespace PhysBAM{

template<class T> class KD_TREE_NODE;
template<class TV> class IMPLICIT_OBJECT;

struct BOX_VISITOR_TRIVIAL
{
    ARRAY<ARRAY<int> >& intersection_list;

    BOX_VISITOR_TRIVIAL(ARRAY<ARRAY<int> >& intersection_list)
        :intersection_list(intersection_list)
    {}

    bool Cull_Self(const int self_box_index) const
    {return false;}

    bool Cull(const int self_box_index,const int other_box_index) const
    {return false;}

    void Store(const int self_box_index,const int other_box_index) const
    {intersection_list(self_box_index).Append(other_box_index);}
};

template<class T_NESTED_VISITOR>
struct BOX_VISITOR_MPI
{
    T_NESTED_VISITOR nested_visitor;
    const ARRAY<char> &processors1,&processors2; // entries are 2*has_ours + has_greater 

    BOX_VISITOR_MPI(T_NESTED_VISITOR nested_visitor,const ARRAY<char>& processors1,const ARRAY<char>& processors2)
        :nested_visitor(nested_visitor),processors1(processors1),processors2(processors2)
    {}

    bool Cull_Self(const int box) const
    {return processors1(box)<2 || nested_visitor.Cull_Self(box);} // has_ours

    bool Cull(const int box1,const int box2) const
    {return (processors1(box1)<<processors2(box2))<4 || nested_visitor.Cull(box1,box2);} // has_ours_1 and (has_ours_2 or has_greater_2) or has_ours_2 and has_greater_1

    void Store(const int self_box_index,const int other_box_index)
    {nested_visitor.Store(self_box_index,other_box_index);}
};

template<class TV>
class BOX_HIERARCHY:public NONCOPYABLE
{
private:
    typedef typename TV::SCALAR T;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::HYPERPLANE T_HYPERPLANE;
public:
    int leaves,root;
    ARRAY<int> parents;
    ARRAY<VECTOR<int,2> > children;
    ARRAY<RANGE<TV> > box_hierarchy;
    mutable STACK<int> traversal_stack;
    mutable STACK<VECTOR<int,2> > dual_traversal_stack;
    ARRAY<T> box_radius;

    BOX_HIERARCHY();
    virtual ~BOX_HIERARCHY();

    bool Leaf(const int box) const
    {return box<=leaves;}

    void Update_Box_Radii()
    {Update_Leaf_Box_Radii();Update_Nonleaf_Box_Radii();}

    void Update_Leaf_Box_Radii()
    {Calculate_Bounding_Box_Radii(box_hierarchy,box_radius);}

    virtual void Intersection_List(const TV& point,ARRAY<int>& intersection_list,const T thickness_over_two=0) const
    {Intersection_List(root,point,intersection_list,thickness_over_two);}
    
    virtual void Intersection_List(const RANGE<TV>& test_box,ARRAY<int>& intersection_list,const T thickness_over_two=0) const
    {Intersection_List(root,test_box,intersection_list,thickness_over_two);}
    
    virtual void Intersection_List(const ORIENTED_BOX<TV>& test_box,ARRAY<int>& intersection_list) const
    {Intersection_List(root,test_box,intersection_list);}

    virtual void Intersection_List(const T_HYPERPLANE& test_plane,ARRAY<int>& intersection_list,const T thickness_over_two=0) const
    {Intersection_List(root,test_plane,intersection_list,thickness_over_two);}

    virtual void Intersection_List(const IMPLICIT_OBJECT<TV>& implicit_object,const MATRIX<T,TV::dimension>& rotation,const TV& translation,ARRAY<int>& intersection_list,const T contour_value=0) const
    {Intersection_List(root,implicit_object,rotation,translation,intersection_list,contour_value);}
    
//#####################################################################
    virtual void Clean_Memory();
protected:
    int Initialize_Hierarchy_Using_KD_Tree_Helper(KD_TREE_NODE<T>* node);
public:
    virtual void Initialize_Hierarchy_Using_KD_Tree();
    void Set_Leaf_Boxes(const ARRAY<RANGE<TV> >& boxes,const bool reinitialize=false);
    void Thicken_Leaf_Boxes(const T extra_thickness);
    void Update_Nonleaf_Boxes();
    void Update_Modified_Nonleaf_Boxes(ARRAY<bool>& modified);
    virtual void Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius);
    void Update_Nonleaf_Box_Radii();
    // for internal use - but octrees use them as well so they're not private
    template<class T_THICKNESS> void Intersection_List(const int box,const RANGE<TV>& test_box,ARRAY<int>& intersection_list,const T_THICKNESS thickness_over_two) const;
protected:
    template<class T_THICKNESS> void Intersection_List(const int box,const TV& point,ARRAY<int>& intersection_list,const T_THICKNESS thickness_over_two) const;
    void Intersection_List(const int box,const ORIENTED_BOX<TV>& test_box,ARRAY<int>& intersection_list) const;
    void Intersection_List(const int box,const T_HYPERPLANE& test_plane,ARRAY<int>& intersection_list,const T thickness_over_two) const;
    void Intersection_List(const int box,const IMPLICIT_OBJECT<TV>& implicit_object,const MATRIX<T,TV::dimension>& rotation,const TV& translation,ARRAY<int>& intersection_list,
        const T contour_value) const;
    template<class T_VISITOR,class T_THICKNESS> void Intersection_List(const BOX_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const int self_box,const int other_box,
        const T_THICKNESS extra_thickness=ZERO()) const;
public:
    template<class T_VISITOR,class T_THICKNESS> void Intersection_List(const BOX_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const T_THICKNESS extra_thickness=ZERO()) const;
    template<class T_VISITOR> void Intersection_List(T_VISITOR& visitor) const;
//#####################################################################
};   
}
#endif

