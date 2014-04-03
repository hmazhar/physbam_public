//#####################################################################
// Copyright 2007, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/IMPLICIT_OBJECT_INTERSECTOR.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_INTERSECTOR<TV>::
IMPLICIT_OBJECT_INTERSECTOR(const IMPLICIT_OBJECT<TV>& implicit_object_input)
    :implicit_object(implicit_object_input),maximum_refinement_depth(0),minimum_refinement_depth(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_INTERSECTOR<TV>::
~IMPLICIT_OBJECT_INTERSECTOR()
{}
//#####################################################################
// Function Negative_Material
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INTERSECTOR<TV>::
Negative_Material_In_Cell(const GRID<TV>& grid,const TV_INT& cell_index,const bool force_full_refinement)
{
    // refinement step
    static TV_INT phi_indices[GRID<TV>::number_of_nodes_per_cell];grid.Nodes_In_Cell_From_Minimum_Corner_Node(cell_index,phi_indices);

    cell_particle_X.Resize(GRID<TV>::number_of_nodes_per_cell);
    for(int i=1;i<=GRID<TV>::number_of_nodes_per_cell;i++) cell_particle_X(i)=grid.Node(phi_indices[i-1]);
    cell_refinement_simplices.Remove_All();Refined_Object_Initialization_Helper(cell_refinement_simplices);

    int last_node=cell_particle_X.m;
    cell_phis.Resize(last_node);
    // compute phis for extant nodes
    T minimum_phi=FLT_MAX,maximum_phi=-FLT_MAX;

    for(int i=1;i<=last_node;i++){T& phi=cell_phis(i);
        phi=grid_nodal_phis(phi_indices[i-1]);
        if(phi<minimum_phi) minimum_phi=phi;
        if(phi>maximum_phi) maximum_phi=phi;}

    if(minimum_phi*maximum_phi>0 && min(abs(minimum_phi),abs(maximum_phi))>grid.Minimum_Edge_Length()) return minimum_phi<=0?full_cell_size:0;

    int unrefined_point_count=last_node;
    int last_parent_simplex=0,first_parent_simplex=1;
    for(int depth=0;depth<maximum_refinement_depth;depth++){
        first_parent_simplex=last_parent_simplex+1;
        last_parent_simplex=cell_refinement_simplices.m;
        for(int s=last_parent_simplex;s>=first_parent_simplex;s--){
            const T_ELEMENT simplex=cell_refinement_simplices(s);
            if(depth < minimum_refinement_depth || Refinement_Condition(cell_particle_X,cell_phis,simplex) || force_full_refinement){
                last_parent_simplex--;
                cell_refinement_simplices.Remove_Index_Lazy(s);
                Refine_Simplex(*this,cell_refinement_simplices,cell_particle_X,simplex);}}
        cell_phis.Resize(cell_particle_X.m);
        // compute phis for extant nodes
        for(int i=last_node+1;i<=cell_phis.m;i++) cell_phis(i)=Phi(cell_particle_X(i));
        last_node=cell_phis.m;}

    if(maximum_refinement_depth>0 && cell_phis.m==unrefined_point_count) return minimum_phi<=0?full_cell_size:0;
    // compute material
    T negative_material=0;
    for(int i=1;i<=cell_refinement_simplices.m;i++) negative_material+=T_SIMPLEX::Negative_Material(cell_particle_X,cell_phis,cell_refinement_simplices(i));
    return clamp(negative_material,(T)0,full_cell_size);
}
//#####################################################################
// Function Negative_Material
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INTERSECTOR<TV>::
Negative_Material_In_Box(const RANGE<TV>& box,const bool force_full_refinement)
{
    T_ARRAYS_VECTOR corners;box.Corners(corners);
    cell_particle_X.Resize(1<<TV::dimension); // redundant
    Make_List(cell_particle_X,corners);
    cell_refinement_simplices.Remove_All();Refined_Object_Initialization_Helper(cell_refinement_simplices);
    T minimum_edge_length=box.Edge_Lengths().Min();

    int last_node=cell_particle_X.m;
    cell_phis.Resize(last_node);
    // compute phis for extant nodes
    T minimum_phi=FLT_MAX,maximum_phi=-FLT_MAX;

    for(int i=1;i<=last_node;i++){T& phi=cell_phis(i);
        phi=Extended_Phi(cell_particle_X(i));
        if(phi<minimum_phi) minimum_phi=phi;
        if(phi>maximum_phi) maximum_phi=phi;}

    if(minimum_phi*maximum_phi>0 && min(abs(minimum_phi),abs(maximum_phi))>minimum_edge_length) return minimum_phi<=0?box.Size():0;

    int unrefined_point_count=last_node;
    int last_parent_simplex=0,first_parent_simplex=1;
    for(int depth=0;depth<maximum_refinement_depth;depth++){
        first_parent_simplex=last_parent_simplex+1;
        last_parent_simplex=cell_refinement_simplices.m;
        for(int s=last_parent_simplex;s>=first_parent_simplex;s--){
            const T_ELEMENT simplex=cell_refinement_simplices(s);
            if(depth < minimum_refinement_depth || Refinement_Condition(cell_particle_X,cell_phis,simplex) || force_full_refinement){
                last_parent_simplex--;
                cell_refinement_simplices.Remove_Index_Lazy(s);
                Refine_Simplex(*this,cell_refinement_simplices,cell_particle_X,simplex);}}
        cell_phis.Resize(cell_particle_X.m);
        // compute phis for extant nodes
        for(int i=last_node+1;i<=cell_phis.m;i++) cell_phis(i)=Extended_Phi(cell_particle_X(i));
        last_node=cell_phis.m;}

    if(maximum_refinement_depth>0 && cell_phis.m==unrefined_point_count) return minimum_phi<=0?box.Size():0;
    // compute material
    T negative_material=0;
    for(int i=1;i<=cell_refinement_simplices.m;i++) negative_material+=T_SIMPLEX::Negative_Material(cell_particle_X,cell_phis,cell_refinement_simplices(i));
    return clamp(negative_material,(T)0,box.Size());
}
//#####################################################################
// Function Negative_Material
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_INTERSECTOR<TV>::
Negative_Material_In_Box_Excluding_Object(const RANGE<TV>& box,const ARRAY<IMPLICIT_OBJECT<TV>*>& excluded_implicit_objects,const bool force_full_refinement)
{
    T_ARRAYS_VECTOR corners;box.Corners(corners);
    cell_particle_X.Resize(1<<TV::dimension); // redundant
    Make_List(cell_particle_X,corners);
    cell_refinement_simplices.Remove_All();Refined_Object_Initialization_Helper(cell_refinement_simplices);
    T minimum_edge_length=box.Edge_Lengths().Min();

    int last_node=cell_particle_X.m;
    cell_phis.Resize(last_node);
    // compute phis for extant nodes
    T minimum_phi=FLT_MAX,maximum_phi=-FLT_MAX;

    for(int i=1;i<=last_node;i++){T& phi=cell_phis(i);
        phi=Extended_Phi(cell_particle_X(i));
        // TODO This is the slow way to do it...
        for(int j=1;j<=excluded_implicit_objects.m;++j) phi=max(phi,-excluded_implicit_objects(j)->Extended_Phi(cell_particle_X(i)));
        if(phi<minimum_phi) minimum_phi=phi;
        if(phi>maximum_phi) maximum_phi=phi;}

    if(minimum_phi*maximum_phi>0 && min(abs(minimum_phi),abs(maximum_phi))>minimum_edge_length) return minimum_phi<=0?box.Size():0;

    int unrefined_point_count=last_node;
    int last_parent_simplex=0,first_parent_simplex=1;
    for(int depth=0;depth<maximum_refinement_depth;depth++){
        first_parent_simplex=last_parent_simplex+1;
        last_parent_simplex=cell_refinement_simplices.m;
        for(int s=last_parent_simplex;s>=first_parent_simplex;s--){
            const T_ELEMENT simplex=cell_refinement_simplices(s);
            if(depth < minimum_refinement_depth || Refinement_Condition(cell_particle_X,cell_phis,simplex) || force_full_refinement){
                last_parent_simplex--;
                cell_refinement_simplices.Remove_Index_Lazy(s);
                Refine_Simplex(*this,cell_refinement_simplices,cell_particle_X,simplex);}}
        cell_phis.Resize(cell_particle_X.m);
        // compute phis for extant nodes
        for(int i=last_node+1;i<=cell_phis.m;i++){T& phi=cell_phis(i);
            phi=Extended_Phi(cell_particle_X(i));
            for(int j=1;j<=excluded_implicit_objects.m;++j) phi=max(phi,-excluded_implicit_objects(j)->Extended_Phi(cell_particle_X(i)));}
        last_node=cell_phis.m;}

    if(maximum_refinement_depth>0 && cell_phis.m==unrefined_point_count) return minimum_phi<=0?box.Size():0;
    // compute material
    T negative_material=0;
    for(int i=1;i<=cell_refinement_simplices.m;i++) negative_material+=T_SIMPLEX::Negative_Material(cell_particle_X,cell_phis,cell_refinement_simplices(i));
    return clamp(negative_material,(T)0,box.Size());
}
//#####################################################################
// Function Iterative_Find_Interface
//#####################################################################
template<class TV> TV IMPLICIT_OBJECT_INTERSECTOR<TV>::
Iterative_Find_Interface(TV left,TV right,const int iterations) const
{
    T phi_left=Extended_Phi(left),phi_right=Extended_Phi(right);
    TV theta=LINEAR_INTERPOLATION<T,TV>::Linear(left,right,LEVELSET_UTILITIES<T>::Theta(phi_left,phi_right));
    int phi_left_sign=(phi_left<=0?-1:1),phi_right_sign=(phi_right<=0?-1:1);
    for(int i=1;i<=iterations;i++){
        T phi=Extended_Phi(theta);
        int phi_sign=(phi<=0?-1:1);
        if(phi_left_sign*phi_sign<0){
            right=theta;phi_right=phi;phi_right_sign=phi_sign;
            theta=LINEAR_INTERPOLATION<T,TV>::Linear(left,theta,LEVELSET_UTILITIES<T>::Theta(phi_left,phi));}
        else if(phi_right_sign*phi_sign<0){
            left=theta;phi_left=phi;phi_left_sign=phi_sign;
            theta=LINEAR_INTERPOLATION<T,TV>::Linear(theta,right,LEVELSET_UTILITIES<T>::Theta(phi,phi_right));}
        else break;}
    return theta;
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template IMPLICIT_OBJECT_INTERSECTOR<VECTOR<T,d> >::IMPLICIT_OBJECT_INTERSECTOR(const IMPLICIT_OBJECT<VECTOR<T,d> >& implicit_object_input); \
    template IMPLICIT_OBJECT_INTERSECTOR<VECTOR<T,d> >::~IMPLICIT_OBJECT_INTERSECTOR(); \
    template T IMPLICIT_OBJECT_INTERSECTOR<VECTOR<T,d> >::Negative_Material_In_Box(const RANGE<VECTOR<T,d> >& box,const bool force_full_refinement); \
    template T IMPLICIT_OBJECT_INTERSECTOR<VECTOR<T,d> >::Negative_Material_In_Box_Excluding_Object(const RANGE<VECTOR<T,d> >& box,const ARRAY<IMPLICIT_OBJECT<VECTOR<T,d> >*>& excluded_implicit_object,const bool force_full_refinement); \
    template VECTOR<T,d> IMPLICIT_OBJECT_INTERSECTOR<VECTOR<T,d> >::Iterative_Find_Interface(VECTOR<T,d> left,VECTOR<T,d> right,const int iterations) const;
INSTANTIATION_HELPER(float,1)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1)
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
#endif
