//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Neil Molino, Andrew Selle, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT
//#####################################################################
#ifndef __IMPLICIT_OBJECT__
#define __IMPLICIT_OBJECT__

#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
namespace PhysBAM{

template<class TV>
class IMPLICIT_OBJECT:public STRUCTURE<TV>
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
    typedef VECTOR<T,d-1> T_CURVATURES;
public:
    typedef TV VECTOR_T;

    BOX<TV> box; // box containing the voxelized implicit surface
    bool use_secondary_interpolation;

    IMPLICIT_OBJECT();
    virtual ~IMPLICIT_OBJECT();

    void Use_Secondary_Interpolation(const bool use_secondary_interpolation_input=true)
    {use_secondary_interpolation=use_secondary_interpolation_input;}

    virtual RANGE<TV>& Box()
    {return box;}

//#####################################################################
    virtual void Update_Box();
    virtual void Update_Minimum_Cell_Size(const int maximum_depth=0);
    virtual T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const;
    virtual T operator()(const TV& location) const;
    virtual T Signed_Distance(const TV& location) const; // to make this class compatible with the geometry classes
    virtual T Extended_Phi(const TV& location) const;
    virtual T Phi_Secondary(const TV& location) const;
    virtual TV Normal(const TV& location,const int aggregate=-1) const;
    virtual TV Extended_Normal(const TV& location,const int aggregate=-1) const;
    virtual void Compute_Normals();
    virtual void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true);
    void Rescale(const T scaling_factor) PHYSBAM_OVERRIDE;
    virtual void Translate(const TV& translation);
    virtual void Inflate(const T inflation_distance);
    virtual bool Inside(const TV& location,const T thickness_over_two=0) const;
    virtual bool Outside(const TV& location,const T thickness_over_two=0) const;
    virtual bool Boundary(const TV& location,const T thickness_over_two) const;
    virtual bool Lazy_Inside(const TV& location,const T contour_value=0) const;
    virtual bool Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value=0) const;
    virtual bool Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const;
    virtual bool Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const;
    virtual bool Lazy_Outside(const TV& location,const T contour_value=0) const;
    virtual bool Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const;
    virtual bool Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const;
    virtual T Min_Phi() const;
    virtual bool Intersection(RAY<TV>& ray,const T thickness=0) const;
    virtual TV Closest_Point_On_Boundary(const TV& location,const T tolerance=0,const int max_iterations=1,T* distance=0) const;
    virtual TV Velocity(const TV& location) const;
    // the following only exist in 3d
    virtual T_SYMMETRIC_MATRIX Hessian(const TV& X) const;
    virtual T_CURVATURES Principal_Curvatures(const TV& X) const;
    virtual T Integration_Step(const T phi) const;
    virtual T Minimum_Cell_Size() const;
//#####################################################################
};
}
#endif
