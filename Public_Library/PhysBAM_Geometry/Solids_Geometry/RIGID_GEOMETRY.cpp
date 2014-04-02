//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Ronald Fedkiw, Jon Gretarsson, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Frank Losasso, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/PROJECTED_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#endif
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_POINT_SIMPLICES_1D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_SEGMENTED_CURVE_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_GEOMETRY<TV>::
RIGID_GEOMETRY(RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,bool create_collision_geometry,int index)
    :implicit_object(0),simplicial_object(0),rigid_geometry_collection(rigid_geometry_collection_input),is_static(false),bounding_box_up_to_date(false),moving_simplex_hierarchy(0),impulse_accumulator(0)
{
    if(index) particle_index=index;
    else particle_index=rigid_geometry_collection.particles.array_collection->Add_Element();
    assert(!rigid_geometry_collection.particles.rigid_geometry(particle_index));
    rigid_geometry_collection.particles.rigid_geometry(particle_index)=dynamic_cast<RIGID_GEOMETRY<TV>*>(this);

    Set_Surface_Roughness();
    Set_Coefficient_Of_Friction();
    if(create_collision_geometry && rigid_geometry_collection.collision_body_list)
        rigid_geometry_collection.collision_body_list->Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(*this),particle_index,true);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_GEOMETRY<TV>::
~RIGID_GEOMETRY()
{
    rigid_geometry_collection.collision_body_list->Remove_Body(rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(particle_index));
    rigid_geometry_collection.particles.structure_ids(particle_index)=VECTOR<int,3>();
    rigid_geometry_collection.particles.rigid_geometry(particle_index)=0;
    delete moving_simplex_hierarchy;
    delete implicit_object;
}
//#####################################################################
// Function Pointwise_Object_Velocity
//#####################################################################
template<class TV> TV RIGID_GEOMETRY<TV>::
Pointwise_Object_Velocity(const TV& X) const
{
    return Pointwise_Object_Velocity(Twist(),this->X(),X);
}
//#####################################################################
// Function Pointwise_Object_Velocity_At_Particle
//#####################################################################
template<class TV> TV RIGID_GEOMETRY<TV>::
Pointwise_Object_Velocity_At_Particle(const TV& X,const int particle_index) const
{
    return Pointwise_Object_Velocity(X);
}
//#####################################################################
// Function Implicit_Geometry_Extended_Value
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_GEOMETRY<TV>::
Implicit_Geometry_Extended_Value(const TV& location) const
{
    return implicit_object->Extended_Phi(location);
}
//#####################################################################
// Function Implicit_Geometry_Normal
//#####################################################################
template<class TV> TV RIGID_GEOMETRY<TV>::
Implicit_Geometry_Normal(const TV& location,const int aggregate) const
{
    return implicit_object->Normal(location,aggregate);
}
//#####################################################################
// Function Implicit_Geometry_Normal
//#####################################################################
template<class TV> TV RIGID_GEOMETRY<TV>::
Implicit_Geometry_Normal(const TV& location,T& phi_value,const int aggregate,const int location_particle_index) const
{
    phi_value=(*implicit_object)(location);
    return implicit_object->Normal(location,aggregate);
}
//#####################################################################
// Function Implicit_Geometry_Lazy_Inside
//#####################################################################
template<class TV> bool RIGID_GEOMETRY<TV>::
Implicit_Geometry_Lazy_Inside(const TV& location,T contour_value) const
{
    return implicit_object->Lazy_Inside(location,contour_value);
}
//#####################################################################
// Function Implicit_Geometry_Lazy_Inside_And_Value
//#####################################################################
template<class TV> bool RIGID_GEOMETRY<TV>::
Implicit_Geometry_Lazy_Inside_And_Value(const TV& location,T& phi,T contour_value) const
{
    return implicit_object->Lazy_Inside_And_Value(location,phi,contour_value);
}
//#####################################################################
// Function Implicit_Geometry_Extended_Normal
//#####################################################################
template<class TV> TV RIGID_GEOMETRY<TV>::
Implicit_Geometry_Extended_Normal(const TV& location,T& phi_value,const int aggregate,const int location_particle_index) const
{
    phi_value=implicit_object->Extended_Phi(location);
    return implicit_object->Extended_Normal(location,aggregate);
}
//#####################################################################
// Function Simplex_Intersection
//#####################################################################
template<class TV> bool RIGID_GEOMETRY<TV>::
Simplex_Intersection(RAY<TV>& ray,const T collision_thickness) const
{
    RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(INTERSECTION::Intersects(object_space_ray,*simplicial_object,collision_thickness)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    return false;
}
//#####################################################################
// Function Add_Structure_Helper
//#####################################################################
template<class T> void
Add_Structure_Helper(RIGID_GEOMETRY<VECTOR<T,1> >& self,STRUCTURE<VECTOR<T,1> >& structure)
{
    typedef VECTOR<T,1> TV;
    if(POINT_SIMPLICES_1D<T>* point_simplices=dynamic_cast<POINT_SIMPLICES_1D<T>*>(&structure)){ // set up acceleration stuctures too
        if(self.simplicial_object) self.Remove_Structure(self.simplicial_object);
        self.simplicial_object=point_simplices;
        if(!self.simplicial_object->bounding_box) self.simplicial_object->Update_Bounding_Box();}
    else if(IMPLICIT_OBJECT<TV>* implicit_object_input=dynamic_cast<IMPLICIT_OBJECT<TV>*>(&structure)){
        if(self.implicit_object) self.Remove_Structure(self.implicit_object->object_space_implicit_object);
        self.implicit_object=new IMPLICIT_OBJECT_TRANSFORMED<TV,RIGID_GEOMETRY<TV> >(implicit_object_input,false,&self);
        self.implicit_object->Update_Box();self.implicit_object->Compute_Cell_Minimum_And_Maximum(false);} // don't recompute cell min/max if already computed
    self.structures.Append(&structure);
}
template<class T> void
Add_Structure_Helper(RIGID_GEOMETRY<VECTOR<T,2> >& self,STRUCTURE<VECTOR<T,2> >& structure)
{
    typedef VECTOR<T,2> TV;
    if(SEGMENTED_CURVE_2D<T>* segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>*>(&structure)){ // set up acceleration stuctures too
        if(self.simplicial_object) self.Remove_Structure(self.simplicial_object);
        self.simplicial_object=segmented_curve;
        if(!self.simplicial_object->bounding_box) self.simplicial_object->Update_Bounding_Box();}
    else if(IMPLICIT_OBJECT<TV>* implicit_object_input=dynamic_cast<IMPLICIT_OBJECT<TV>*>(&structure)){
        if(self.implicit_object) self.Remove_Structure(self.implicit_object->object_space_implicit_object);
        self.implicit_object=new IMPLICIT_OBJECT_TRANSFORMED<TV,RIGID_GEOMETRY<TV> >(implicit_object_input,false,&self);
        self.implicit_object->Update_Box();self.implicit_object->Compute_Cell_Minimum_And_Maximum(false);} // don't recompute cell min/max if already computed
    else if(TRIANGULATED_AREA<T>* triangulated_area=dynamic_cast<TRIANGULATED_AREA<T>*>(&structure)){
        if(TRIANGULATED_AREA<T>* old_triangulated_area=self.template Find_Structure<TRIANGULATED_AREA<T>*>()) self.Remove_Structure(old_triangulated_area);
        if(!triangulated_area->bounding_box) triangulated_area->Update_Bounding_Box();}
    self.structures.Append(&structure);
}
template<class T> void
Add_Structure_Helper(RIGID_GEOMETRY<VECTOR<T,3> >& self,STRUCTURE<VECTOR<T,3> >& structure)
{
    typedef VECTOR<T,3> TV;
    if(TRIANGULATED_SURFACE<T>* triangulated_surface_input=dynamic_cast<TRIANGULATED_SURFACE<T>*>(&structure)){ // set up acceleration stuctures too
        if(self.simplicial_object) self.Remove_Structure(self.simplicial_object);
        self.simplicial_object=triangulated_surface_input;
        if(!self.simplicial_object->hierarchy){self.simplicial_object->Initialize_Hierarchy();self.simplicial_object->hierarchy->Update_Box_Radii();}
        if(!self.simplicial_object->bounding_box) self.simplicial_object->Update_Bounding_Box();}
    else if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* implicit_object_input=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(&structure)){
        if(self.implicit_object) self.Remove_Structure(self.implicit_object->object_space_implicit_object);
        self.implicit_object=new IMPLICIT_OBJECT_TRANSFORMED<TV,RIGID_GEOMETRY<TV> >(implicit_object_input,false,&self);
        self.implicit_object->Update_Box();self.implicit_object->Compute_Cell_Minimum_And_Maximum(false);} // don't recompute cell min/max if already computed
    else if(IMPLICIT_OBJECT<TV>* implicit_object_input=dynamic_cast<IMPLICIT_OBJECT<TV>*>(&structure)){
        if(self.implicit_object) self.Remove_Structure(self.implicit_object->object_space_implicit_object);
        self.implicit_object=new IMPLICIT_OBJECT_TRANSFORMED<TV,RIGID_GEOMETRY<TV> >(implicit_object_input,false,&self);
        self.implicit_object->Update_Box();self.implicit_object->Compute_Cell_Minimum_And_Maximum(false);} // don't recompute cell min/max if already computed
    else if(TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(&structure)){
        if(TETRAHEDRALIZED_VOLUME<T>* old_tetrahedralized_volume=self.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>()) self.Remove_Structure(old_tetrahedralized_volume);
        if(!tetrahedralized_volume->bounding_box) tetrahedralized_volume->Update_Bounding_Box();}
    self.structures.Append(&structure);
}
//#####################################################################
// Function Add_Structure
//#####################################################################
template<class TV> void RIGID_GEOMETRY<TV>::
Add_Structure(STRUCTURE<TV>& structure)
{
    Add_Structure_Helper(dynamic_cast<RIGID_GEOMETRY<TV>&>(*this),structure);
}
//#####################################################################
// Function Print_Names
//#####################################################################
template<class TV> void RIGID_GEOMETRY<TV>::
Print_Names(const RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection,const ARRAY<int>& ids)
{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    LOG::cout<<"{";
    for(int i=1;i<=ids.m;i++){LOG::cout<<"\""<<rigid_geometry_collection.particles.rigid_geometry(ids(i))->name<<"\"";if(i<ids.m) LOG::cout<<", ";}
    LOG::cout<<"}";
#endif
}
//#####################################################################
// Function Object_Space_Bounding_Box
//#####################################################################
template<class TV> const RANGE<TV>& RIGID_GEOMETRY<TV>::
Object_Space_Bounding_Box() const // at least one of triangulated surface or implicit surface must exist for this to work
{
    if(simplicial_object){PHYSBAM_ASSERT(simplicial_object->bounding_box);return *simplicial_object->bounding_box;}
    else{PHYSBAM_ASSERT(implicit_object);return implicit_object->object_space_implicit_object->Box();}
}
//#####################################################################
// Function Remove_Structure
//#####################################################################
template<class TV> void RIGID_GEOMETRY<TV>::
Remove_Structure(STRUCTURE<TV>* structure)
{
    structures.Remove_Index_Lazy(structures.Find(structure));
}
//#####################################################################
// Function World_Space_Simplex_Bounding_Box
//#####################################################################
template<class TV> RANGE<TV> RIGID_GEOMETRY<TV>::
World_Space_Simplex_Bounding_Box(const int id) const
{
    const VECTOR<int,TV::dimension>& elements=simplicial_object->mesh.elements(id);
    VECTOR<TV,TV::dimension> pts;
    for(int i=1;i<=TV::dimension;i++) pts(i)=World_Space_Point(simplicial_object->particles.X(elements(i)));
    return RANGE<TV>::Bounding_Box(pts);
}
//#####################################################################
// Function World_Space_Simplex
//#####################################################################
template<class TV> typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX RIGID_GEOMETRY<TV>::
World_Space_Simplex(const int id) const
{
    const VECTOR<int,TV::dimension>& elements=simplicial_object->mesh.elements(id);
    VECTOR<TV,TV::dimension> pts;
    for(int i=1;i<=TV::dimension;i++) pts(i)=World_Space_Point(simplicial_object->particles.X(elements(i)));
    return typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX(pts);
}
//#####################################################################
// Function World_Space_Simplex_Bounding_Box
//#####################################################################
template<class TV> RANGE<TV> RIGID_GEOMETRY<TV>::
World_Space_Simplex_Bounding_Box(const int id,const FRAME<TV>& frame) const
{
    const VECTOR<int,TV::dimension>& elements=simplicial_object->mesh.elements(id);
    VECTOR<TV,TV::dimension> pts;
    for(int i=1;i<=TV::dimension;i++) pts(i)=frame*simplicial_object->particles.X(elements(i));
    return RANGE<TV>::Bounding_Box(pts);
}
template<class TV> void RIGID_GEOMETRY<TV>::
Update_Bounding_Box()
{
    bounding_box_up_to_date=true;
    const RANGE<TV>& box=Object_Space_Bounding_Box();
    oriented_box=T_ORIENTED_BOX(box,Frame());
    axis_aligned_bounding_box=oriented_box.Axis_Aligned_Bounding_Box();
}
template<class TV> void RIGID_GEOMETRY<TV>::
Update_Bounding_Box_From_Implicit_Geometry()
{
    assert(implicit_object);
    bounding_box_up_to_date=true;
    oriented_box=T_ORIENTED_BOX(implicit_object->object_space_implicit_object->box,Frame());
    axis_aligned_bounding_box=oriented_box.Axis_Aligned_Bounding_Box();
}
//#####################################################################
// Function Compute_Velocity_Between_States
//#####################################################################
template<class TV> void RIGID_GEOMETRY<TV>::
Compute_Velocity_Between_States(const RIGID_GEOMETRY_STATE<TV>& state1,const RIGID_GEOMETRY_STATE<TV>& state2,RIGID_GEOMETRY_STATE<TV>& result_state)
{
    RIGID_GEOMETRY_STATE<TV>::Compute_Velocity_Between_States(state1,state2,result_state);
}
//#####################################################################
// Function Interpolate_Between_States
//#####################################################################
template<class TV> void RIGID_GEOMETRY<TV>::
Interpolate_Between_States(const RIGID_GEOMETRY_STATE<TV>& state1,const RIGID_GEOMETRY_STATE<TV>& state2,const T time,RIGID_GEOMETRY_STATE<TV>& interpolated_state)
{
    PHYSBAM_ASSERT(((time>=state1.time && time<=state2.time) || (time<=state1.time && time>=state2.time)) && state1.time<state2.time);
    T alpha=(time-state1.time)/(state2.time-state1.time);alpha=clamp(alpha,(T)0,(T)1);
    interpolated_state.frame.t=TV::Interpolate(state1.frame.t,state2.frame.t,alpha);
    interpolated_state.frame.r=ROTATION<TV>::Spherical_Linear_Interpolation(state1.frame.r,state2.frame.r,alpha);
    interpolated_state.twist.linear=(1-alpha)*state1.twist.linear+alpha*state2.twist.linear;
    interpolated_state.time=time;
}
//#####################################################################
template class RIGID_GEOMETRY<VECTOR<float,1> >;
template class RIGID_GEOMETRY<VECTOR<float,2> >;
template class RIGID_GEOMETRY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_GEOMETRY<VECTOR<double,1> >;
template class RIGID_GEOMETRY<VECTOR<double,2> >;
template class RIGID_GEOMETRY<VECTOR<double,3> >;
#endif
}
