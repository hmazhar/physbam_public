//#####################################################################
// Copyright 2002-2007, Doug Enright, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Neil Molino, Andrew Selle, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#ifndef __LEVELSET_IMPLICIT_OBJECT__
#define __LEVELSET_IMPLICIT_OBJECT__

#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class LEVELSET_IMPLICIT_OBJECT:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
    enum WORKAROUND {d=TV::m};
    typedef VECTOR<T,d-1> T_PRINCIPAL_CURVATURES;
public:
    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;

    T_LEVELSET levelset;
    T_ARRAYS_VECTOR* V;
    INTERPOLATION_UNIFORM<GRID<TV>,TV>* velocity_interpolation;
private:
    static LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,TV> default_velocity_interpolation;
protected:
    T minimum_cell_size;
    bool need_destroy_data;
public:

    LEVELSET_IMPLICIT_OBJECT(GRID<TV>& grid_input,T_ARRAYS_SCALAR& phi_input);
    virtual ~LEVELSET_IMPLICIT_OBJECT();

//###########################################################################
    static LEVELSET_IMPLICIT_OBJECT<TV>* Create();
    void Set_Custom_Interpolation(INTERPOLATION_UNIFORM<GRID<TV>,T>& interpolation);
    void Set_Custom_Secondary_Interpolation(INTERPOLATION_UNIFORM<GRID<TV>,T>& interpolation);
    void Set_Custom_Normal_Interpolation(INTERPOLATION_UNIFORM<GRID<TV>,TV>& interpolation);
    void Set_Custom_Velocity_Interpolation(INTERPOLATION_UNIFORM<GRID<TV>,TV>& interpolation);
    void Update_Box() PHYSBAM_OVERRIDE;
    void Update_Minimum_Cell_Size(const int maximum_depth=0) PHYSBAM_OVERRIDE;
    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const PHYSBAM_OVERRIDE;
    T operator()(const TV& location) const PHYSBAM_OVERRIDE;
    T Extended_Phi(const TV& location) const PHYSBAM_OVERRIDE;
    T Phi_Secondary(const TV& location) const PHYSBAM_OVERRIDE;
    TV Normal(const TV& location,const int aggregate=-1) const PHYSBAM_OVERRIDE;
    TV Extended_Normal(const TV& location,const int aggregate=-1) const PHYSBAM_OVERRIDE;
    void Compute_Normals() PHYSBAM_OVERRIDE;
    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true) PHYSBAM_OVERRIDE;
    void Inflate(const T inflation_distance) PHYSBAM_OVERRIDE;
    bool Lazy_Inside(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Lazy_Outside(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE;
    T Min_Phi() const PHYSBAM_OVERRIDE;
    TV Velocity(const TV& location) const PHYSBAM_OVERRIDE;
    T_SYMMETRIC_MATRIX Hessian(const TV& X) const PHYSBAM_OVERRIDE;
    VECTOR<T,d-1> Principal_Curvatures(const TV& X) const PHYSBAM_OVERRIDE;
    void Rescale(const T scaling_factor) PHYSBAM_OVERRIDE;
    void Translate(const TV& translation) PHYSBAM_OVERRIDE;
    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    virtual std::string Extension() const PHYSBAM_OVERRIDE {return Static_Extension();}
    static std::string Static_Name();
    static std::string Static_Extension();
    T Integration_Step(const T phi) const PHYSBAM_OVERRIDE;
    T Minimum_Cell_Size() const PHYSBAM_OVERRIDE;
//###########################################################################
};
}
#endif
