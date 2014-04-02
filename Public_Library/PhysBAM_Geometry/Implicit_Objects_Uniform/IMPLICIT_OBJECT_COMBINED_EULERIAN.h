//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_COMBINED_EULERIAN
//#####################################################################
#ifndef __IMPLICIT_OBJECT_COMBINED_EULERIAN__
#define __IMPLICIT_OBJECT_COMBINED_EULERIAN__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>

namespace PhysBAM{

template<class TV>
class IMPLICIT_OBJECT_COMBINED_EULERIAN:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;
    typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
    enum WORKAROUND {d=TV::m};
public:
    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;

    IMPLICIT_OBJECT<TV> *implicit_object1,*implicit_object2;
    bool owns_implicit_object1,owns_implicit_object2;
    T alpha;T dt; // alpha is the time step fraction, dt is the time step from implicit_object1 to implicit_object2, i.e. the time you want to interpolate at is dt*alpha
    
    IMPLICIT_OBJECT_COMBINED_EULERIAN(IMPLICIT_OBJECT<TV>* implicit_object1_input,bool owns_implicit_object1_input,IMPLICIT_OBJECT<TV>* implicit_object2_input, bool owns_implicit_object2_input)
        :implicit_object1(implicit_object1_input),implicit_object2(implicit_object2_input),owns_implicit_object1(owns_implicit_object1_input),owns_implicit_object2(owns_implicit_object2_input)
    {}

    virtual ~IMPLICIT_OBJECT_COMBINED_EULERIAN()
    {if(owns_implicit_object1) delete implicit_object1;if(owns_implicit_object2) delete implicit_object2;}

    void Set_Weights(T alpha_input)
    {alpha=alpha_input;}

    void Update_Box() PHYSBAM_OVERRIDE
    {implicit_object1->Update_Box();implicit_object2->Update_Box();box=implicit_object1->box;
    box.Enlarge_Nonempty_Box_To_Include_Points(implicit_object2->box.min_corner,implicit_object2->box.max_corner);}

    void Update_Minimum_Cell_Size(const int maximum_depth=0) PHYSBAM_OVERRIDE
    {implicit_object1->Update_Minimum_Cell_Size(maximum_depth);implicit_object2->Update_Minimum_Cell_Size(maximum_depth);}

    // TODO: box is not used.  Is it still needed?
    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const PHYSBAM_OVERRIDE
    {return min(implicit_object1->Minimum_Cell_Size_Within_Box(box),implicit_object2->Minimum_Cell_Size_Within_Box(box));} 

    TV Velocity(const TV& location) const PHYSBAM_OVERRIDE
    {return (1-alpha)*implicit_object1->Velocity(location)+alpha*implicit_object2->Velocity(location);}

    T operator()(const TV& location) const PHYSBAM_OVERRIDE
    {return Value(location);}

    T Value(const TV& location) const
    {TV V_half=Velocity(location); // half velocity for semi lagrangian
    T phi1=(*implicit_object1)(location-alpha*dt*V_half),phi2=(*implicit_object2)(location+dt*(1-alpha)*V_half); 
    return (1-alpha)*phi1+alpha*phi2;}

    T Extended_Value(const TV& location) const
    {TV V_half=Velocity(location); // half velocity for semi lagrangian
    T phi1=implicit_object1->Extended_Phi(location-alpha*dt*V_half),phi2=implicit_object2->Extended_Phi(location+dt*(1-alpha)*V_half);
    return (1-alpha)*phi1+alpha*phi2;}

    T Signed_Distance(const TV& location) const PHYSBAM_OVERRIDE
    {return (*this)(location);} // to make this class compatible with the geometry classes

    T Extended_Phi(const TV& location) const PHYSBAM_OVERRIDE
    {return Extended_Value(location);}

    T Phi_Secondary(const TV& location) const PHYSBAM_OVERRIDE
    {TV V_half=Velocity(location); // half velocity for semi lagrangian
    T phi1=implicit_object1->Phi_Secondary(location-alpha*dt*V_half),phi2=implicit_object2->Phi_Secondary(location+dt*(1-alpha)*V_half); 
    return (1-alpha)*phi1+alpha*phi2;}

    TV Normal(const TV& location,const int aggregate) const PHYSBAM_OVERRIDE
    {T dx=Minimum_Cell_Size();TV normal;T one_over_2_dx=(T)1/((T)2*dx);
    for(int i=1;i<=TV::m;i++){
        TV offset=dx*TV::Axis_Vector(i);
        normal[i]=(Value(location+offset)-Value(location-offset))*one_over_2_dx;}
    return normal.Normalized();}

    TV Extended_Normal(const TV& location,const int aggregate=-1) const PHYSBAM_OVERRIDE
    {T dx=Minimum_Cell_Size();TV normal;T one_over_2_dx=(T)1/((T)2*dx);
    for(int i=1;i<=TV::m;i++){
        TV offset=dx*TV::Axis_Vector(i);
        normal[i]=(Extended_Value(location+offset)-Extended_Value(location-offset))*one_over_2_dx;}
    return normal.Normalized();}

    void Compute_Normals() PHYSBAM_OVERRIDE
    {} 

    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true) PHYSBAM_OVERRIDE
    {implicit_object1->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);implicit_object2->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);}
    
    void Rescale(const T scaling_factor) PHYSBAM_OVERRIDE
    {implicit_object1->Rescale(scaling_factor);implicit_object2->Rescale(scaling_factor);}

    void Inflate(const T inflation_distance) PHYSBAM_OVERRIDE
    {implicit_object1->Inflate(inflation_distance);implicit_object2->Inflate(inflation_distance);}

    bool Lazy_Inside(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return box.Inside(location,-contour_value) && (*this)(location)<=contour_value;}

    bool Lazy_Inside_And_Value(const TV& location,T& phi,const T contour_value=0) const PHYSBAM_OVERRIDE
    {if(box.Inside(location,-contour_value)){phi=(*this)(location);if(phi<=contour_value) return true;} return false;}

    bool Lazy_Inside_Extended_Levelset(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return box.Inside(location,-contour_value) && Extended_Value(location)<=contour_value;}

    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE
    {if(box.Inside(location,-contour_value)){phi_value=Extended_Phi(location);if(phi_value<=contour_value) return true;} return false;}

    bool Lazy_Outside(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return !Lazy_Inside(location,contour_value);}

    bool Lazy_Outside_Extended_Levelset(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return !Lazy_Inside_Extended_Levelset(location,contour_value);}

    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return !Lazy_Inside_Extended_Levelset_And_Value(location,phi_value,contour_value);}

    T Min_Phi() const PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}

    T_SYMMETRIC_MATRIX Hessian(const TV& X) const PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}

    VECTOR<T,d-1> Principal_Curvatures(const TV& X) const PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}

    T Integration_Step(const T phi) const PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}

    T Minimum_Cell_Size() const PHYSBAM_OVERRIDE
    {return min(implicit_object1->Minimum_Cell_Size(),implicit_object2->Minimum_Cell_Size());}

    bool Intersection(RAY<TV>& ray,const T thickness) const PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################
};
}
#endif

