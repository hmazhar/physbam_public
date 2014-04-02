//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_TRANSFORMED
//#####################################################################
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace PhysBAM{
bool Register_Implicit_Object_Transformed(){
    static bool done=false;if(done) return true;done=true;
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,2>,FRAME<VECTOR<float,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,3>,FRAME<VECTOR<float,3> > > >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,2>,FRAME<VECTOR<double,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,3>,FRAME<VECTOR<double,3> > > >();
#endif
    return true;
}
bool registered_asdlfkjasdlfkjas=Register_Implicit_Object_Transformed();
}
using namespace PhysBAM;
template<class TV,class TRANSFORM> IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<TV>* object_space_implicit_object_input,bool owns_implicit_object_input,const TRANSFORM* transform_input)
    :IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM>(transform_input),owns_implicit_object(owns_implicit_object_input),object_space_implicit_object(object_space_implicit_object_input)
{
}
template<class TV,class TRANSFORM> IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<TV>* object_space_implicit_object_input,bool owns_implicit_object_input,const TV& center_input,T scale_input)
    :IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM>(center_input,scale_input),owns_implicit_object(owns_implicit_object_input),object_space_implicit_object(object_space_implicit_object_input)
{
}
template<class TV,class TRANSFORM> IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
~IMPLICIT_OBJECT_TRANSFORMED()
{
    if(owns_implicit_object) delete object_space_implicit_object;
}
template<class TV,class TRANSFORM> RANGE<TV>& IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Box()
{
    Update_Box();
    return box;
}
template<class TV,class TRANSFORM> void IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Update_Box()
{
    object_space_implicit_object->Update_Box();
    box=World_Space_Box(object_space_implicit_object->box);
}
template<class TV,class TRANSFORM> void IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Update_Minimum_Cell_Size(const int maximum_depth)
{
    object_space_implicit_object->Update_Minimum_Cell_Size(maximum_depth);
}
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const
{
    return World_Space_Length(object_space_implicit_object->Minimum_Cell_Size_Within_Box(Object_Space_Box(object_space_implicit_object->box)));
}
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
operator()(const TV& location) const
{
    return World_Space_Length((*object_space_implicit_object)(Object_Space_Point(location)));
}
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Signed_Distance(const TV& location) const
{
    return (*this)(location); // to make this class compatible with the geometry classes
}
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Extended_Phi(const TV& location) const
{
    return Extended_Value(location);
}
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Phi_Secondary(const TV& location) const
{
    return World_Space_Length(object_space_implicit_object->Phi_Secondary(Object_Space_Point(location)));
}
template<class TV,class TRANSFORM> TV IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Normal(const TV& location,const int aggregate) const
{
    return World_Space_Unitless_Vector(object_space_implicit_object->Normal(Object_Space_Point(location),aggregate));
}
template<class TV,class TRANSFORM> TV IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Extended_Normal(const TV& location,const int aggregate) const
{
    return World_Space_Unitless_Vector(object_space_implicit_object->Extended_Normal(Object_Space_Point(location),aggregate));
}
template<class TV,class TRANSFORM> void IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Compute_Normals()
{
    object_space_implicit_object->Compute_Normals();
}
template<class TV,class TRANSFORM> void IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    object_space_implicit_object->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);
}
template<class TV,class TRANSFORM> void IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Rescale(const T scaling_factor)
{
    if(!Scale_Transform(scaling_factor)) object_space_implicit_object->Rescale(scaling_factor);
}
template<class TV,class TRANSFORM> void IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Inflate(const T inflation_distance)
{
    object_space_implicit_object->Inflate(Object_Space_Length(inflation_distance));
}
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Inside(const TV& location,const T thickness_over_two) const
{
    return object_space_implicit_object->Inside(Object_Space_Point(location),Object_Space_Length(thickness_over_two));
}
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Inside(const TV& location,const T contour_value) const
{
    return object_space_implicit_object->Lazy_Inside(Object_Space_Point(location),Object_Space_Length(contour_value));
}
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Inside_And_Value(const TV& location,T& phi,const T contour_value) const
{
    bool result=object_space_implicit_object->Lazy_Inside_And_Value(Object_Space_Point(location),phi,Object_Space_Length(contour_value));
    phi=World_Space_Length(phi);
    return result;
}
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Inside_Extended_Levelset(const TV& location,const T contour_value) const
{
    return object_space_implicit_object->Lazy_Inside_Extended_Levelset(Object_Space_Point(location),Object_Space_Length(contour_value));
}
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Inside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value) const
{
    bool result=object_space_implicit_object->Lazy_Inside_Extended_Levelset_And_Value(Object_Space_Point(location),phi_value,Object_Space_Length(contour_value));
    phi_value=World_Space_Length(phi_value);
    return result;
}
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Outside(const TV& location,const T thickness_over_two) const
{
    return object_space_implicit_object->Outside(Object_Space_Point(location),Object_Space_Length(thickness_over_two));
}
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Outside(const TV& location,const T contour_value) const
{
    return object_space_implicit_object->Lazy_Outside(Object_Space_Point(location),Object_Space_Length(contour_value));
}
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Outside_Extended_Levelset(const TV& location,const T contour_value) const
{
    return object_space_implicit_object->Lazy_Outside_Extended_Levelset(Object_Space_Point(location),Object_Space_Length(contour_value));
}
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Lazy_Outside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value) const
{
    bool result=object_space_implicit_object->Lazy_Outside_Extended_Levelset_And_Value(Object_Space_Point(location),phi_value,Object_Space_Length(contour_value));
    phi_value=World_Space_Length(phi_value);
    return result;
}
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Min_Phi() const
{
    return World_Space_Length(object_space_implicit_object->Min_Phi());
}
template<class TV,class TRANSFORM> TV IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Velocity(const TV& location) const
{
    return World_Space_Vector(object_space_implicit_object->Velocity(Object_Space_Point(location)));
}
template<class TV,class TRANSFORM> typename IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::T_SYMMETRIC_MATRIX IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Hessian(const TV& X) const
{
    return World_Space_Length_Hessian(object_space_implicit_object->Hessian(Object_Space_Point(X)));
}
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Integration_Step(const T phi) const
{
    return World_Space_Length(object_space_implicit_object->Integration_Step(Object_Space_Length(phi)));
}
template<class TV,class TRANSFORM> typename TV::SCALAR IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Minimum_Cell_Size() const
{
    return World_Space_Length(object_space_implicit_object->Minimum_Cell_Size());
}
template<class TV,class TRANSFORM> bool IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Intersection(RAY<TV>& ray,const T thickness) const
{
    RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(object_space_implicit_object->Intersection(object_space_ray,thickness)){
        ray.semi_infinite=false;
        ray.t_max=World_Space_Length(object_space_ray.t_max);
        ray.aggregate_id=object_space_ray.aggregate_id;
        return true;}
    else return false;
}
template<class TV,class TRANSFORM> TV IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Closest_Point_On_Boundary(const TV& location,const T tolerance,const int max_iterations,T* distance) const
{
    TV result=World_Space_Point(object_space_implicit_object->Closest_Point_On_Boundary(Object_Space_Point(location),Object_Space_Length(tolerance),max_iterations,distance));
    if(distance) *distance=World_Space_Length(*distance);
    return result;
}
template<class TV,class TRANSFORM> std::string IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::
Name() const
{
    return Static_Name();
}
//#####################################################################
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,1>,RIGID_GEOMETRY<VECTOR<float,1> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,1> >*,bool,RIGID_GEOMETRY<VECTOR<float,1> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,2>,RIGID_GEOMETRY<VECTOR<float,2> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,2> >*,bool,RIGID_GEOMETRY<VECTOR<float,2> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,3>,RIGID_GEOMETRY<VECTOR<float,3> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,3> >*,bool,RIGID_GEOMETRY<VECTOR<float,3> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,1>,FRAME<VECTOR<float,1> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,1> >*,bool,FRAME<VECTOR<float,1> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,2>,FRAME<VECTOR<float,2> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,2> >*,bool,FRAME<VECTOR<float,2> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,3>,FRAME<VECTOR<float,3> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,3> >*,bool,FRAME<VECTOR<float,3> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,1>,float>::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,1> >*,bool,VECTOR<float,1> const&,const float);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,2>,float>::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,2> >*,bool,VECTOR<float,2> const&,const float);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<float,3>,float>::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<float,3> >*,bool,VECTOR<float,3> const&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,1>,RIGID_GEOMETRY<VECTOR<double,1> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,1> >*,bool,RIGID_GEOMETRY<VECTOR<double,1> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,2>,RIGID_GEOMETRY<VECTOR<double,2> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,2> >*,bool,RIGID_GEOMETRY<VECTOR<double,2> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,3>,RIGID_GEOMETRY<VECTOR<double,3> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,3> >*,bool,RIGID_GEOMETRY<VECTOR<double,3> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,1>,FRAME<VECTOR<double,1> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,1> >*,bool,FRAME<VECTOR<double,1> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,2>,FRAME<VECTOR<double,2> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,2> >*,bool,FRAME<VECTOR<double,2> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,3>,FRAME<VECTOR<double,3> > >::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,3> >*,bool,FRAME<VECTOR<double,3> > const*);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,1>,double>::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,1> >*,bool,VECTOR<double,1> const&,const double);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,2>,double>::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,2> >*,bool,VECTOR<double,2> const&,const double);
template IMPLICIT_OBJECT_TRANSFORMED<VECTOR<double,3>,double>::IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<VECTOR<double,3> >*,bool,VECTOR<double,3> const&,const double);
#endif
