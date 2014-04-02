//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Avi Robinson-Mosher, Craig Schroeder, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ORIENTED_BOX
//#####################################################################
#ifndef __ORIENTED_BOX__
#define __ORIENTED_BOX__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class TV>
class ORIENTED_BOX
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND{d=TV::dimension};
    typedef typename MATRIX_POLICY<TV>::DIAGONAL_MATRIX T_DIAGONAL_MATRIX;
public:
    typedef TV VECTOR_T;

    TV corner; // root corner of the box
    MATRIX<T,d> edges; // principal edges of the box, emanating from the corner

    ORIENTED_BOX()
        :corner(TV()),edges(MATRIX<T,d>::Identity_Matrix())
    {}

    ORIENTED_BOX(const TV& corner_input,const TV& edge1_input,const TV& edge2_input,const TV& edge3_input)
        :corner(corner_input),edges(edge1_input,edge2_input,edge3_input)
    {STATIC_ASSERT(d==3);}

    ORIENTED_BOX(const TV& corner_input,const TV& edge1_input,const TV& edge2_input)
        :corner(corner_input),edges(edge1_input,edge2_input)
    {STATIC_ASSERT(d==2);}

    ORIENTED_BOX(const TV& corner_input,const MATRIX<T,d>& edges_input)
        :corner(corner_input),edges(edges_input)
    {}

    ORIENTED_BOX(const RANGE<TV>& box,const ROTATION<TV>& rotation)
    {
        corner=rotation.Rotate(box.min_corner);
        edges=rotation.Rotation_Matrix()*T_DIAGONAL_MATRIX(box.Edge_Lengths());
    }

    ORIENTED_BOX(const RANGE<TV>& box,const ROTATION<TV>& rotation,const TV& corner_input)
        :corner(corner_input)
    {
        edges=rotation.Rotation_Matrix()*T_DIAGONAL_MATRIX(box.Edge_Lengths());
    }

    ORIENTED_BOX(const RANGE<TV>& box,const FRAME<TV>& frame)
    {
        corner=frame*box.min_corner;
        edges=frame.r.Rotation_Matrix()*T_DIAGONAL_MATRIX(box.Edge_Lengths());
    }

    ORIENTED_BOX(const RANGE<TV>& box,const MATRIX<T,d+1>& transform)
    {
        corner=transform.Homogeneous_Times(box.Minimum_Corner());
        edges=transform.Extract_Rotation()*T_DIAGONAL_MATRIX(box.Edge_Lengths());
    }

    TV Center() const
    {return corner+(T).5*edges.Column_Sum();}

    void Scale_About_Center(const T factor)
    {TV center=Center(); // compute before modifying edges
    edges*=factor;corner=center-(T).5*edges.Column_Sum();}

    void Scale_About_Center(const TV factor)
    {TV center=Center(); // compute before modifying edges
    edges*=T_DIAGONAL_MATRIX(factor);corner=center-(T).5*edges.Column_Sum();}

    ORIENTED_BOX Scaled_About_Center(const T factor) const
    {ORIENTED_BOX box(*this);box.Scale_About_Center(factor);return box;}

    ORIENTED_BOX Scaled_About_Center(const TV edge_factor)
    {ORIENTED_BOX box(*this);box.Scale_About_Center(edge_factor);return box;}

    RANGE<TV> Axis_Aligned_Bounding_Box() const
    {RANGE<TV> box(corner);for(int i=1;i<=d;i++) box.Enlarge_By_Sign(edges.Column(i));return box;}

    RANGE<TV> Bounding_Box() const // for templatization purposes
    {return Axis_Aligned_Bounding_Box();}

    bool Lazy_Inside(const TV &location) const
    {TV vec=location-corner,edge_projection=edges.Transpose_Times(vec);
    // TODO: you presently cannot be inside a degenerate box
    for(int i=1;i<=d;i++) if(edge_projection(i)<0 || edge_projection(i)>edges.Column(i).Magnitude_Squared() || edges.Column(i).Magnitude_Squared()==0) return false;
    return true;}

    void Generate_Unit_Direction(VECTOR<T,1>& v,int i) const
    {v=VECTOR<T,1>(edges(1,1)<0?-1:1);}

    void Generate_Unit_Direction(VECTOR<T,2>& v,int i) const
    {v=edges.Column(i);if(!v.Normalize()) v=edges.Column(3-i).Perpendicular().Normalized();}

    void Generate_Unit_Direction(VECTOR<T,3>& v,int i) const
    {v=edges.Column(i);if(!v.Normalize()) v=VECTOR<T,3>::Cross_Product(edges.Column(i==1?2:1),edges.Column(i==3?2:3)).Normalized();}

    // Handles a zero-length edge
    ORIENTED_BOX Thickened(const T thickness) const
    {ORIENTED_BOX box(*this);TV u;
    for(int i=1;i<=d;i++){
        Generate_Unit_Direction(u,i);
        u*=thickness;
        box.corner-=u;
        box.edges.Column(i)+=u+u;}
    return box;}

    VECTOR<T,TV::dimension-1> Principal_Curvatures(const TV& X) const
    {return VECTOR<T,TV::dimension-1>();}

//#####################################################################
    bool Intersection(const ORIENTED_BOX& box) const;
    bool Intersection(const RANGE<TV>& box) const;
    T Signed_Distance(const TV& X) const;
    TV Normal(const TV& X) const;
    bool Separating_Test(const ORIENTED_BOX& box,const TV& plane_normal_direction) const;
    bool Separating_Test(const RANGE<TV>& box,const TV& plane_normal_direction) const;
    static std::string Name();
private:
    void Project_Points_Onto_Line(const TV& direction,T& line_min,T& line_max) const;
//#####################################################################
};
}
#endif
