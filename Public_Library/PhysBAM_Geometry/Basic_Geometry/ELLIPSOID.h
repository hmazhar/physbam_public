//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELLIPSOID
//#####################################################################
#ifndef __ELLIPSOID__
#define __ELLIPSOID__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
namespace PhysBAM{

template<class T>
class ELLIPSOID
{
    typedef VECTOR<T,3> TV;
public:
    TV center;
    DIAGONAL_MATRIX<T,3> radii;
    ROTATION<TV> orientation;

    ELLIPSOID()
        :radii(1,1,1)
    {}

    ELLIPSOID(const TV& center_input,const DIAGONAL_MATRIX<T,3>& radii_input,const ROTATION<TV>& orientation_input=ROTATION<TV>())
        :center(center_input),radii(radii_input),orientation(orientation_input.Normalized())
    {}

    ELLIPSOID(const TV& center_input,const DIAGONAL_MATRIX<T,3>& radii_input,const TV& x_axis,const TV& y_axis,const TV& z_axis) // assumes axes are orthogonal
        :center(center_input),radii(radii_input)
    {
        Initialize(x_axis.Normalized(),y_axis.Normalized(),z_axis.Normalized());
    }

    ELLIPSOID(const TV& center_input,TV scaled_x_axis,TV scaled_y_axis,TV scaled_z_axis) // assumes axes are orthogonal
        :center(center_input),radii(scaled_x_axis.Normalize(),scaled_y_axis.Normalize(),scaled_z_axis.Normalize())
    {
        Initialize(scaled_x_axis,scaled_y_axis,scaled_z_axis);
    }

    ELLIPSOID(const TV& center_input,const DIAGONAL_MATRIX<T,3>& radii_input,const MATRIX<T,3>& axes) // assumes axes are orthonormal
        :center(center_input),radii(radii_input)
    {
        Initialize(axes.Column(1).Normalized(),axes.Column(2).Normalized(),axes.Column(3).Normalized());
    }

private:
    void Initialize(const TV& x_axis,const TV& y_axis,const TV& z_axis)
    {orientation=ROTATION<TV>(MATRIX<T,3>::Rotation_Matrix(x_axis,y_axis,z_axis)).Normalized();}
public:

    ORIENTED_BOX<TV> Oriented_Bounding_Box() const
    {MATRIX<T,3> axes(orientation.Rotation_Matrix()*radii);return ORIENTED_BOX<TV>(center-axes.Column_Sum(),(T)2*axes.Column(1),(T)2*axes.Column(2),(T)2*axes.Column(3));}

    RANGE<TV> Bounding_Box() const
    {return Oriented_Bounding_Box().Axis_Aligned_Bounding_Box();}

    SYMMETRIC_MATRIX<T,3> Metric_Tensor() const
    {return (orientation.Rotation_Matrix()*radii.Inverse()).Outer_Product_Matrix();}

    T Lipschitz_Constant() const
    {return 1/radii.Min();}

    T Volume() const
    {return (T)four_thirds_pi*radii.Determinant();}

//#####################################################################
    TV Normal(const TV& location) const;
    bool Inside(const TV& location,const T thickness_over_two=0) const;
    bool Outside(const TV& location,const T thickness_over_two=0) const;
    bool Boundary(const TV& location,const T thickness_over_two) const;
    TV Approximate_Surface(const TV& location) const; // not necessarily the closest point on the surface, but approximates it
    T Approximate_Signed_Distance(const TV& location) const; // not the true signed distance, but has correct inside/outside
    template<class T_ARRAY_TV> static ELLIPSOID<T> Covariance_Ellipsoid(const T_ARRAY_TV& points);
//#####################################################################
};
}
#endif
