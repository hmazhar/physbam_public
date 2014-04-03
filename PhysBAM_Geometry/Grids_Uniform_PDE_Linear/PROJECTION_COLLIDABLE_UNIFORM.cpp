//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/PROJECTION_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PROJECTION_COLLIDABLE_UNIFORM<T_GRID>::
PROJECTION_COLLIDABLE_UNIFORM(const T_GRID& mac_grid,const bool multiphase,const bool use_poisson,const bool use_variable_beta,THREAD_QUEUE* thread_queue)
    :laplace_collidable(0),poisson_collidable(0)
{
    if(use_poisson){
        poisson=poisson_collidable=new POISSON_COLLIDABLE_UNIFORM<T_GRID>(p_grid,p,true,multiphase,true);
        if(use_variable_beta) poisson->Set_Variable_beta();elliptic_solver=poisson;collidable_solver=poisson_collidable;}
    else{
        laplace=laplace_collidable=new LAPLACE_COLLIDABLE_UNIFORM<T_GRID>(p_grid,p,true,false,true,thread_queue);elliptic_solver=laplace;collidable_solver=laplace_collidable;}
    Initialize_Grid(mac_grid);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PROJECTION_COLLIDABLE_UNIFORM<T_GRID>::
PROJECTION_COLLIDABLE_UNIFORM(const T_GRID& mac_grid,T_LEVELSET& levelset_input)
    :laplace_collidable(0),poisson_collidable(0)
{
    poisson=poisson_collidable=new POISSON_COLLIDABLE_UNIFORM<T_GRID>(p_grid,p,levelset_input,true,false,true);elliptic_solver=poisson;collidable_solver=poisson_collidable;
    Initialize_Grid(mac_grid);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> PROJECTION_COLLIDABLE_UNIFORM<T_GRID>::
~PROJECTION_COLLIDABLE_UNIFORM()
{
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class T_GRID> void PROJECTION_COLLIDABLE_UNIFORM<T_GRID>::
Initialize_Grid(const T_GRID& mac_grid)
{
    BASE::Initialize_Grid(mac_grid);
}
//#####################################################################
// Function Apply_Pressure
//#####################################################################
template<class T_GRID> void PROJECTION_COLLIDABLE_UNIFORM<T_GRID>::
Apply_Pressure(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,bool scale_by_dt)
{
    // find divergence free u, v and w 
    if(collidable_solver->second_order_cut_cell_method){
        //Zero_Out_Neumann_Pocket_Velocities(face_velocities); TODO: Why is this here?
        T_ARRAYS_BOOL& psi_D=elliptic_solver->psi_D;
        T_FACE_ARRAYS_BOOL& psi_N=elliptic_solver->psi_N;
        T_ARRAYS_SCALAR& phi=collidable_solver->levelset->phi;
        TV dx=p_grid.dX,one_over_dx=Inverse(dx);
        if(scale_by_dt) p*=dt;
        if(laplace){
            for(FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
                if(!psi_N.Component(axis)(face_index) && !(psi_D(first_cell) && psi_D(second_cell))){
                    if(psi_D(first_cell) && !psi_D(second_cell) && LEVELSET_UTILITIES<T>::Interface(phi(second_cell),phi(first_cell)))
                        face_velocities.Component(axis)(face_index)-=(p(second_cell)-collidable_solver->u_interface.Component(axis)(face_index))
                            /(max(LEVELSET_UTILITIES<T>::Theta(phi(second_cell),phi(first_cell)),collidable_solver->second_order_cut_cell_threshold)*dx[axis]);
                    else if(!psi_D(first_cell) && psi_D(second_cell) && LEVELSET_UTILITIES<T>::Interface(phi(first_cell),phi(second_cell)))
                        face_velocities.Component(axis)(face_index)-=(collidable_solver->u_interface.Component(axis)(face_index)-p(first_cell))
                            /(max(LEVELSET_UTILITIES<T>::Theta(phi(first_cell),phi(second_cell)),collidable_solver->second_order_cut_cell_threshold)*dx[axis]);
                    else face_velocities.Component(axis)(face_index)-=(p(second_cell)-p(first_cell))*one_over_dx[axis];}}}
        else if(poisson){
            assert(!poisson->use_variable_beta); // assumes constant beta in each phase
            for(FACE_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
                if(!psi_N.Component(axis)(face_index) && !(psi_D(first_cell) && psi_D(second_cell))){
                    if(psi_D(first_cell) && !psi_D(second_cell) && LEVELSET_UTILITIES<T>::Interface(phi(second_cell),phi(first_cell)))
                        face_velocities.Component(axis)(face_index)-=poisson->beta_face.Component(axis)(face_index)*(p(second_cell)-collidable_solver->u_interface.Component(axis)(face_index))
                            /(max(LEVELSET_UTILITIES<T>::Theta(phi(second_cell),phi(first_cell)),collidable_solver->second_order_cut_cell_threshold)*dx[axis]);
                    else if(!psi_D(first_cell) && psi_D(second_cell) && LEVELSET_UTILITIES<T>::Interface(phi(first_cell),phi(second_cell)))
                        face_velocities.Component(axis)(face_index)-=poisson->beta_face.Component(axis)(face_index)*(collidable_solver->u_interface.Component(axis)(face_index)-p(first_cell))
                            /(max(LEVELSET_UTILITIES<T>::Theta(phi(first_cell),phi(second_cell)),collidable_solver->second_order_cut_cell_threshold)*dx[axis]);
                    else face_velocities.Component(axis)(face_index)-=poisson->beta_face.Component(axis)(face_index)*(p(second_cell)-p(first_cell))*one_over_dx[axis];}}}
            if(scale_by_dt) p*=1/dt;}
    else
        BASE::Apply_Pressure(face_velocities,dt,time,scale_by_dt);
}
//#####################################################################
template class PROJECTION_COLLIDABLE_UNIFORM<GRID<VECTOR<float,1> > >;
template class PROJECTION_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> > >;
template class PROJECTION_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PROJECTION_COLLIDABLE_UNIFORM<GRID<VECTOR<double,1> > >;
template class PROJECTION_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> > >;
template class PROJECTION_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
