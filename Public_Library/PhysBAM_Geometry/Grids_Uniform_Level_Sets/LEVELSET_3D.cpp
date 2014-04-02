//#####################################################################
// Copyright 2002-2008, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Frank Losasso, Neil Molino, Igor Neverov, Avi Robinson-Mosher,
//     Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LEVELSET_3D<T_GRID>::
LEVELSET_3D(T_GRID& grid_input,ARRAY<T,TV_INT>& phi_input,const int number_of_ghost_cells_input)
    :LEVELSET_UNIFORM<T_GRID>(grid_input,phi_input,number_of_ghost_cells_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> LEVELSET_3D<T_GRID>::
~LEVELSET_3D()
{}
//#####################################################################
// Function Hessian
//#####################################################################
template<class T_GRID> SYMMETRIC_MATRIX<typename T_GRID::SCALAR,3> LEVELSET_3D<T_GRID>::
Hessian(const TV& X) const
{
    T one_over_dx=1/grid.dX.x,one_over_dy=1/grid.dX.y,one_over_dz=1/grid.dX.z;T two_phi_center=2*Phi(X);
    T phi_xx=(Phi(TV(X.x+grid.dX.x,X.y,X.z))-two_phi_center+Phi(TV(X.x-grid.dX.x,X.y,X.z)))*sqr(one_over_dx),
       phi_yy=(Phi(TV(X.x,X.y+grid.dX.y,X.z))-two_phi_center+Phi(TV(X.x,X.y-grid.dX.y,X.z)))*sqr(one_over_dy),
       phi_zz=(Phi(TV(X.x,X.y,X.z+grid.dX.z))-two_phi_center+Phi(TV(X.x,X.y,X.z-grid.dX.z)))*sqr(one_over_dz),
       phi_xy=(Phi(TV(X.x+grid.dX.x,X.y+grid.dX.y,X.z))-Phi(TV(X.x+grid.dX.x,X.y-grid.dX.y,X.z))
                  -Phi(TV(X.x-grid.dX.x,X.y+grid.dX.y,X.z))+Phi(TV(X.x-grid.dX.x,X.y-grid.dX.y,X.z)))*(T).25*one_over_dx*one_over_dy,
       phi_xz=(Phi(TV(X.x+grid.dX.x,X.y,X.z+grid.dX.z))-Phi(TV(X.x+grid.dX.x,X.y,X.z-grid.dX.z))
                  -Phi(TV(X.x-grid.dX.x,X.y,X.z+grid.dX.z))+Phi(TV(X.x-grid.dX.x,X.y,X.z-grid.dX.z)))*(T).25*one_over_dx*one_over_dz,
       phi_yz=(Phi(TV(X.x,X.y+grid.dX.y,X.z+grid.dX.z))-Phi(TV(X.x,X.y+grid.dX.y,X.z-grid.dX.z))
                   -Phi(TV(X.x,X.y-grid.dX.y,X.z+grid.dX.z))+Phi(TV(X.x,X.y-grid.dX.y,X.z-grid.dX.z)))*(T).25*one_over_dy*one_over_dz;
    return SYMMETRIC_MATRIX<T,3>(phi_xx,phi_xy,phi_xz,phi_yy,phi_yz,phi_zz);
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T_GRID> VECTOR<typename T_GRID::SCALAR,2> LEVELSET_3D<T_GRID>::
Principal_Curvatures(const TV& X) const
{
    TV grad_phi=TV((Phi(TV(X.x+grid.dX.x,X.y,X.z))-Phi(TV(X.x-grid.dX.x,X.y,X.z)))/(2*grid.dX.x),
                                     (Phi(TV(X.x,X.y+grid.dX.y,X.z))-Phi(TV(X.x,X.y-grid.dX.y,X.z)))/(2*grid.dX.y),
                                     (Phi(TV(X.x,X.y,X.z+grid.dX.z))-Phi(TV(X.x,X.y,X.z-grid.dX.z)))/(2*grid.dX.z));
    TV N=grad_phi;T grad_phi_magnitude=N.Normalize();
    SYMMETRIC_MATRIX<T,3> P=(T)1-SYMMETRIC_MATRIX<T,3>::Outer_Product(N),M=SYMMETRIC_MATRIX<T,3>::Conjugate(P,Hessian(X))/grad_phi_magnitude;
    T trace=M.Trace();
    QUADRATIC<T> quadratic(-1,trace,sqr(M(2,1))-M(1,1)*M(2,2)+sqr(M(3,1))-M(1,1)*M(3,3)+sqr(M(3,2))-M(2,2)*M(3,3));
    quadratic.Compute_Roots();
    if(quadratic.roots == 0) (T).5*VECTOR<T,2>(trace,trace);
    else if(quadratic.roots == 1) return VECTOR<T,2>(quadratic.root1,quadratic.root1);
    return VECTOR<T,2>(quadratic.root1,quadratic.root2);
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
// note that sqrt(phix^2+phiy^2+phiz^2)=1 if it's a distance function
template<class T_GRID> void LEVELSET_3D<T_GRID>::
Compute_Normals(const T time)
{
    int ghost_cells=3;
    T one_over_two_dx=2*grid.one_over_dX.x,one_over_two_dy=2*grid.one_over_dX.y,one_over_two_dz=2*grid.one_over_dX.z;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    int z=phi_ghost.counts.z,yz=phi_ghost.counts.y*z;

    if(!normals) normals=new ARRAY<TV,VECTOR<int,3> >(grid.Domain_Indices(ghost_cells-1));
    for(CELL_ITERATOR iterator(grid,ghost_cells-1);iterator.Valid();iterator.Next()){
        const TV_INT& cell = iterator.Cell_Index();
        int index=phi_ghost.Standard_Index(iterator.Cell_Index());
        VECTOR<T,3>& N((*normals)(cell));
        N.x=(phi_ghost.array(index+yz)-phi_ghost.array(index-yz))*one_over_two_dx;
        N.y=(phi_ghost.array(index+z)-phi_ghost.array(index-z))*one_over_two_dy;
        N.z=(phi_ghost.array(index+1)-phi_ghost.array(index-1))*one_over_two_dz;
        N.Normalize();}
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
// kappa = - DIV(normal), negative for negative phi inside, positive for positive phi inside, sqrt(phix^2+phiy^2+phiy^2)=1 for distance functions
template<class T_GRID> void LEVELSET_3D<T_GRID>::
Compute_Curvature(const T time)
{
    int ghost_cells=3;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);

    if(!curvature) curvature=new ARRAY<T,TV_INT>(grid.Domain_Indices(ghost_cells-1));
    for(CELL_ITERATOR iterator(grid,ghost_cells-1);iterator.Valid();iterator.Next())
        (*curvature)(iterator.Cell_Index())=Compute_Curvature(phi_ghost,iterator.Cell_Index());
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_3D<T_GRID>::
Compute_Curvature(const ARRAY<T,TV_INT>& phi_input,const VECTOR<int,3>& index) const
{
    T one_over_two_dx=1/(2*grid.dX.x),one_over_two_dy=1/(2*grid.dX.y),one_over_two_dz=1/(2*grid.dX.z);
    T one_over_dx_squared=1/sqr(grid.dX.x),one_over_dy_squared=1/sqr(grid.dX.y),one_over_dz_squared=1/sqr(grid.dX.z);
    T one_over_four_dx_dy=1/(4*grid.dX.x*grid.dX.y),one_over_four_dx_dz=1/(4*grid.dX.x*grid.dX.z),one_over_four_dy_dz=1/(4*grid.dX.y*grid.dX.z);
    T max_curvature=1/grid.min_dX; // max resolution
    int z=phi_input.counts.z,yz=phi_input.counts.y*z;

    int s_index=phi_input.Standard_Index(index);
    T phix=(phi_input.array(s_index+yz)-phi_input.array(s_index-yz))*one_over_two_dx;
    T phixx=(phi_input.array(s_index+yz)-2*phi_input.array(s_index)+phi_input.array(s_index-yz))*one_over_dx_squared;
    T phiy=(phi_input.array(s_index+z)-phi_input.array(s_index-z))*one_over_two_dy;
    T phiyy=(phi_input.array(s_index+z)-2*phi_input.array(s_index)+phi_input.array(s_index-z))*one_over_dy_squared;
    T phiz=(phi_input.array(s_index+1)-phi_input.array(s_index-1))*one_over_two_dz;
    T phizz=(phi_input.array(s_index+1)-2*phi_input.array(s_index)+phi_input.array(s_index-1))*one_over_dz_squared;
    T phixy=(phi_input.array(s_index+yz+z)-phi_input.array(s_index+yz-z)-phi_input.array(s_index-yz+z)+phi_input.array(s_index-yz-z))*one_over_four_dx_dy;
    T phixz=(phi_input.array(s_index+yz+1)-phi_input.array(s_index+yz-1)-phi_input.array(s_index-yz+1)+phi_input.array(s_index-yz-1))*one_over_four_dx_dz;
    T phiyz=(phi_input.array(s_index+z+1)-phi_input.array(s_index+z-1)-phi_input.array(s_index-z+1)+phi_input.array(s_index-z-1))*one_over_four_dy_dz;
    T denominator=sqrt(sqr(phix)+sqr(phiy)+sqr(phiz)),curvature;
    if(denominator >= small_number)
        curvature=-(sqr(phix)*phiyy-2*phix*phiy*phixy+sqr(phiy)*phixx+sqr(phix)*phizz-2*phix*phiz*phixz+sqr(phiz)*phixx+sqr(phiy)*phizz-2*phiy*phiz*phiyz+sqr(phiz)*phiyy)/cube(denominator);
    else curvature=LEVELSET_UTILITIES<T>::Sign(phi_input(index))*max_curvature;
    return minmag(curvature,LEVELSET_UTILITIES<T>::Sign(curvature)*max_curvature);
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_3D<T_GRID>::
Compute_Curvature(const TV& location) const
{
    // TODO: optimize
    ARRAY<T,TV_INT> phi_stencil(-1,1,-1,1,-1,1);
    for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++) for(int k=-1;k<=1;k++){
        TV phi_location=grid.dX.x*TV((T)i,0,0)+grid.dX.y*TV(0,(T)j,0)+grid.dX.z*TV(0,0,(T)k)+location;
        phi_stencil(i,j,k)=Phi(phi_location);}
    return Compute_Curvature(phi_stencil,VECTOR<int,3>(0,0,0));
}
//#####################################################################
// Function Compute_Cell_Minimum_And_Maximum
//#####################################################################
template<class T_GRID> void LEVELSET_3D<T_GRID>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    if(!recompute_if_exists && cell_range) return;
    if(!cell_range) cell_range=new ARRAY<RANGE<VECTOR<T,1> >,VECTOR<int,3> >(phi.domain.min_corner.x,phi.domain.max_corner.x-1,phi.domain.min_corner.y,phi.domain.max_corner.y-1,phi.domain.min_corner.z,phi.domain.max_corner.z-1);
    TV_INT i;
    for(i.x=phi.domain.min_corner.x;i.x<=phi.domain.max_corner.x-1;i.x++) for(i.y=phi.domain.min_corner.y;i.y<=phi.domain.max_corner.y-1;i.y++) for(i.z=phi.domain.min_corner.z;i.z<=phi.domain.max_corner.z-1;i.z++){
        int index=phi.Standard_Index(i);
        int z=phi.counts.z,yz=index+phi.counts.y*z;
        T phi1=phi.array(index),phi2=phi.array(yz),phi3=phi.array(index+z),phi4=phi.array(index+1),
           phi5=phi.array(yz+z),phi6=phi.array(yz+1),phi7=phi.array(index+z+1),
           phi8=phi.array(yz+z+1);
        (*cell_range)(i)=RANGE<VECTOR<T,1> >(min(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8),max(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8));}
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T_GRID> void LEVELSET_3D<T_GRID>::
Fast_Marching_Method(const T time,const T stopping_distance,const ARRAY<VECTOR<int,3> >* seed_indices,const bool add_seed_indices_for_ghost_cells)
{
    Get_Signed_Distance_Using_FMM(phi,time,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells);
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T_GRID> void LEVELSET_3D<T_GRID>::
Fast_Marching_Method(ARRAY<bool,VECTOR<int,3> >& seed_indices,const T time,const T stopping_distance,const bool add_seed_indices_for_ghost_cells)
{
    Get_Signed_Distance_Using_FMM(phi,seed_indices,time,stopping_distance,add_seed_indices_for_ghost_cells);
}
//#####################################################################
// Function Get_Signed_Distance_Using_FMM
//#####################################################################
template<class T_GRID> void LEVELSET_3D<T_GRID>::
Get_Signed_Distance_Using_FMM(ARRAY<T,TV_INT>& signed_distance,const T time,const T stopping_distance,const ARRAY<VECTOR<int,3> >* seed_indices,const bool add_seed_indices_for_ghost_cells)
{
    const int ghost_cells=max(2*number_of_ghost_cells+1,1-phi.Domain_Indices().Minimum_Corner()(1));
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells),false);boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    FAST_MARCHING_METHOD_UNIFORM<T_GRID > fmm(*this,ghost_cells,thread_queue);
    fmm.Fast_Marching_Method(phi_ghost,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells);
    ARRAY<T,TV_INT>::Get(signed_distance,phi_ghost);
    boundary->Apply_Boundary_Condition(grid,signed_distance,time);
}
//#####################################################################
// Function Get_Signed_Distance_Using_FMM
//#####################################################################
template<class T_GRID> void LEVELSET_3D<T_GRID>::
Get_Signed_Distance_Using_FMM(ARRAY<T,TV_INT>& signed_distance,ARRAY<bool,VECTOR<int,3> >& seed_indices,const T time,const T stopping_distance,const bool add_seed_indices_for_ghost_cells)
{
    const int ghost_cells=max(2*number_of_ghost_cells+1,1-phi.Domain_Indices().Minimum_Corner()(1));
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells),false);boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    FAST_MARCHING_METHOD_UNIFORM<T_GRID > fmm(*this,ghost_cells,thread_queue);
    fmm.Fast_Marching_Method(phi_ghost,seed_indices,stopping_distance,add_seed_indices_for_ghost_cells);
    ARRAY<T,TV_INT>::Get(signed_distance,phi_ghost);
    boundary->Apply_Boundary_Condition(grid,signed_distance,time);
}
//#####################################################################
// Function Fast_Marching_Method_Outside_Band
//#####################################################################
template<class T_GRID> void LEVELSET_3D<T_GRID>::
Fast_Marching_Method_Outside_Band(const T half_band_width,const T time,const T stopping_distance)
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;
    int ghost_cells=number_of_ghost_cells;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells),false);boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    FAST_MARCHING_METHOD_UNIFORM<T_GRID > fmm(*this,ghost_cells,thread_queue);
    fmm.Fast_Marching_Method(phi_ghost,stopping_distance);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++) if(abs(phi_ghost(i,j,ij)) > half_band_width) phi(i,j,ij)=phi_ghost(i,j,ij);
    boundary->Apply_Boundary_Condition(grid,phi,time);
}
//#####################################################################
// Function Approximate_Surface_Area
//#####################################################################
// calculates the approximate perimeter using delta functions
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_3D<T_GRID>::
Approximate_Surface_Area(const T interface_thickness,const T time) const
{
    int ghost_cells=number_of_ghost_cells;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    T interface_half_width=interface_thickness*grid.dX.Max()/2,one_over_two_dx=1/(2*grid.dX.x),one_over_two_dy=1/(2*grid.dX.y),one_over_two_dz=1/(2*grid.dX.z),surface_area=0;
    for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int ij=1;ij<=grid.counts.z;ij++)
        surface_area+=(LEVELSET_UTILITIES<T>::Delta(phi_ghost(i,j,ij),interface_half_width)*sqrt(sqr((phi_ghost(i+1,j,ij)-phi_ghost(i-1,j,ij))*one_over_two_dx)+
                                sqr((phi_ghost(i,j+1,ij)-phi_ghost(i,j-1,ij))*one_over_two_dy)+sqr((phi_ghost(i,j,ij+1)-phi_ghost(i,j,ij-1))*one_over_two_dz)));
    return surface_area*grid.dX.x*grid.dX.y*grid.dX.z;
}
//#####################################################################
// Function Calculate_Triangulated_Surface_From_Marching_Tetrahedra
//#####################################################################
// uses levelset grid for tet marching - faster than version below because we don't need to interpolate phi
template<class T_GRID> void LEVELSET_3D<T_GRID>::
Calculate_Triangulated_Surface_From_Marching_Tetrahedra(TRIANGULATED_SURFACE<T>& triangulated_surface,const bool include_ghost_values) const
{
    triangulated_surface.Clean_Memory();triangulated_surface.mesh.Clean_Memory();triangulated_surface.particles.array_collection->Clean_Memory();
    int m_start=1,m_end=grid.counts.x,n_start=1,n_end=grid.counts.y,mn_start=1,mn_end=grid.counts.z;
    if(include_ghost_values){m_start=phi.domain.min_corner.x;m_end=phi.domain.max_corner.x;n_start=phi.domain.min_corner.y;n_end=phi.domain.max_corner.y;mn_start=phi.domain.min_corner.z;mn_end=phi.domain.max_corner.z;}
    ARRAY<VECTOR<int,6>,VECTOR<int,3> > edge(m_start,m_end,n_start,n_end,mn_start,mn_end);
    // create particles
    for(int i=m_start;i<=m_end;i++) for(int j=n_start;j<=n_end;j++) for(int k=mn_start;k<=mn_end;k++){TV_INT index(i,j,k);
        if(i<m_end) edge(index)(1)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,index,VECTOR<int,3>(i+1,j,k));
        if(j<n_end) edge(index)(2)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,index,VECTOR<int,3>(i,j+1,k));
        if(k<mn_end)edge(index)(3)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,index,VECTOR<int,3>(i,j,k+1));
        if((i+j+k)%2 == 0){
            if(j<n_end&&k<mn_end)edge(index)(4)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,VECTOR<int,3>(i,j+1,k),VECTOR<int,3>(i,j,k+1));
            if(i<m_end&&k<mn_end)edge(index)(5)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,VECTOR<int,3>(i+1,j,k),VECTOR<int,3>(i,j,k+1));
            if(i<m_end&&j<n_end) edge(index)(6)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,VECTOR<int,3>(i,j+1,k),VECTOR<int,3>(i+1,j,k));}
        else{
            if(j<n_end&&k<mn_end)edge(index)(4)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,index,VECTOR<int,3>(i,j+1,k+1));
            if(i<m_end&&k<mn_end)edge(index)(5)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,index,VECTOR<int,3>(i+1,j,k+1));
            if(i<m_end&&j<n_end) edge(index)(6)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,index,VECTOR<int,3>(i+1,j+1,k));}}
    // calculate triangles
    for(int i=m_start;i<=m_end-1;i++) for(int j=n_start;j<=n_end-1;j++) for(int k=mn_start;k<=mn_end-1;k++)
        if((i+j+k)%2 == 0){
            Append_Triangles(triangulated_surface,edge(i,j,k)(1),edge(i,j,k)(2),edge(i,j,k)(3),edge(i,j,k)(6),edge(i,j,k)(4),edge(i,j,k)(5),phi(i,j,k)); // bottom left
            Append_Triangles(triangulated_surface,edge(i,j,k+1)(6),edge(i+1,j,k)(4),edge(i+1,j,k+1)(2),edge(i,j,k)(5),edge(i+1,j,k)(3),edge(i,j,k+1)(1),phi(i+1,j+1,k+1)); // bottom right
            Append_Triangles(triangulated_surface,edge(i,j+1,k)(1),edge(i+1,j,k)(2),edge(i+1,j+1,k)(3),edge(i,j,k)(6),edge(i+1,j,k)(4),edge(i,j+1,k)(5),phi(i+1,j+1,k)); // top front
            Append_Triangles(triangulated_surface,edge(i,j+1,k)(3),edge(i,j,k)(4),edge(i,j+1,k)(5),edge(i,j,k+1)(2),edge(i,j,k+1)(6),edge(i,j+1,k+1)(1),phi(i,j+1,k)); // top back
            Append_Triangles(triangulated_surface,edge(i,j,k)(6),edge(i,j+1,k)(5),edge(i,j,k)(4),edge(i+1,j,k)(4),edge(i,j,k+1)(6),edge(i,j,k)(5),phi(i,j+1,k));} // center
        else{
            Append_Triangles(triangulated_surface,edge(i,j,k)(6),edge(i+1,j,k)(2),edge(i+1,j,k)(4),edge(i,j,k)(1),edge(i+1,j,k)(3),edge(i,j,k)(5),phi(i+1,j+1,k)); // bottom front
            Append_Triangles(triangulated_surface,edge(i,j,k)(5),edge(i,j,k)(4),edge(i,j,k)(3),edge(i,j,k+1)(6),edge(i,j,k+1)(2),edge(i,j,k+1)(1),phi(i,j,k)); // bottom back
            Append_Triangles(triangulated_surface,edge(i,j,k)(6),edge(i,j,k)(2),edge(i,j,k)(4),edge(i,j+1,k)(1),edge(i,j+1,k)(3),edge(i,j+1,k)(5),phi(i,j,k)); // top left
            Append_Triangles(triangulated_surface,edge(i,j+1,k)(5),edge(i,j+1,k+1)(1),edge(i,j,k+1)(6),edge(i+1,j+1,k)(3),edge(i+1,j,k+1)(2),edge(i+1,j,k)(4),phi(i,j+1,k+1)); // top right
            Append_Triangles(triangulated_surface,edge(i,j,k)(6),edge(i,j,k)(4),edge(i,j,k)(5),edge(i,j+1,k)(5),edge(i,j,k+1)(6),edge(i+1,j,k)(4),phi(i,j,k));} // center
    triangulated_surface.mesh.number_nodes=triangulated_surface.particles.array_collection->Size();
    triangulated_surface.Remove_Degenerate_Triangles();
}
//#####################################################################
// Function If_Zero_Crossing_Add_Particle
//#####################################################################
template<class T_GRID> int LEVELSET_3D<T_GRID>::
If_Zero_Crossing_Add_Particle_By_Index(TRIANGULATED_SURFACE<T>& triangulated_surface,const VECTOR<int,3>& index1,const VECTOR<int,3>& index2) const
{
    int index=0;T phi1=phi(index1),phi2=phi(index2);
    if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
        index=triangulated_surface.particles.array_collection->Add_Element();
        triangulated_surface.particles.X(index)=LINEAR_INTERPOLATION<T,TV>::Linear(grid.X(index1),grid.X(index2),LEVELSET_UTILITIES<T>::Theta(phi1,phi2));}
    return index;
}
//#####################################################################
// Function Calculate_Triangulated_Surface_From_Marching_Tetrahedra
//#####################################################################
template<class T_GRID> void LEVELSET_3D<T_GRID>::
Calculate_Triangulated_Surface_From_Marching_Tetrahedra(const T_GRID& tet_grid,TRIANGULATED_SURFACE<T>& triangulated_surface) const
{
    assert(tet_grid.domain.min_corner.x >= grid.domain.min_corner.x && tet_grid.domain.max_corner.x <= grid.domain.max_corner.x && tet_grid.domain.min_corner.y >= grid.domain.min_corner.y && tet_grid.domain.max_corner.y <= grid.domain.max_corner.y && tet_grid.domain.min_corner.z >= grid.domain.min_corner.z &&
               tet_grid.domain.max_corner.z <= grid.domain.max_corner.z);
    triangulated_surface.Clean_Memory();triangulated_surface.mesh.Clean_Memory();triangulated_surface.particles.array_collection->Clean_Memory();
    ARRAY<VECTOR<int,6>,VECTOR<int,3> > edge(1,tet_grid.counts.x,1,tet_grid.counts.y,1,tet_grid.counts.z);
    // create particles
    int i;for(i=1;i<=tet_grid.counts.x;i++) for(int j=1;j<=tet_grid.counts.y;j++) for(int k=1;k<=tet_grid.counts.z;k++){
        edge(i,j,k)(1)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j,k),tet_grid.X(i+1,j,k));
        edge(i,j,k)(2)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j,k),tet_grid.X(i,j+1,k));
        edge(i,j,k)(3)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j,k),tet_grid.X(i,j,k+1));
        if((i+j+k)%2 == 0){
            edge(i,j,k)(4)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j+1,k),tet_grid.X(i,j,k+1));
            edge(i,j,k)(5)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i+1,j,k),tet_grid.X(i,j,k+1));
            edge(i,j,k)(6)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j+1,k),tet_grid.X(i+1,j,k));}
        else{
            edge(i,j,k)(4)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j,k),tet_grid.X(i,j+1,k+1));
            edge(i,j,k)(5)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j,k),tet_grid.X(i+1,j,k+1));
            edge(i,j,k)(6)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j,k),tet_grid.X(i+1,j+1,k));}}
    // calculate triangles
    for(i=1;i<tet_grid.counts.x;i++) for(int j=1;j<tet_grid.counts.y;j++) for(int k=1;k<tet_grid.counts.z;k++)
        if((i+j+k)%2 == 0){
            Append_Triangles(triangulated_surface,edge(i,j,k)(1),edge(i,j,k)(2),edge(i,j,k)(3),edge(i,j,k)(6),edge(i,j,k)(4),edge(i,j,k)(5),Phi(tet_grid.X(i,j,k))); // bottom left
            Append_Triangles(triangulated_surface,edge(i,j,k+1)(6),edge(i+1,j,k)(4),edge(i+1,j,k+1)(2),edge(i,j,k)(5),edge(i+1,j,k)(3),edge(i,j,k+1)(1),Phi(tet_grid.X(i+1,j+1,k+1))); // bottom right
            Append_Triangles(triangulated_surface,edge(i,j+1,k)(1),edge(i+1,j,k)(2),edge(i+1,j+1,k)(3),edge(i,j,k)(6),edge(i+1,j,k)(4),edge(i,j+1,k)(5),Phi(tet_grid.X(i+1,j+1,k))); // top front
            Append_Triangles(triangulated_surface,edge(i,j+1,k)(3),edge(i,j,k)(4),edge(i,j+1,k)(5),edge(i,j,k+1)(2),edge(i,j,k+1)(6),edge(i,j+1,k+1)(1),Phi(tet_grid.X(i,j+1,k))); // top back
            Append_Triangles(triangulated_surface,edge(i,j,k)(6),edge(i,j+1,k)(5),edge(i,j,k)(4),edge(i+1,j,k)(4),edge(i,j,k+1)(6),edge(i,j,k)(5),Phi(tet_grid.X(i,j+1,k)));} // center
        else{
            Append_Triangles(triangulated_surface,edge(i,j,k)(6),edge(i+1,j,k)(2),edge(i+1,j,k)(4),edge(i,j,k)(1),edge(i+1,j,k)(3),edge(i,j,k)(5),Phi(tet_grid.X(i+1,j+1,k))); // bottom front
            Append_Triangles(triangulated_surface,edge(i,j,k)(5),edge(i,j,k)(4),edge(i,j,k)(3),edge(i,j,k+1)(6),edge(i,j,k+1)(2),edge(i,j,k+1)(1),Phi(tet_grid.X(i,j,k))); // bottom back
            Append_Triangles(triangulated_surface,edge(i,j,k)(6),edge(i,j,k)(2),edge(i,j,k)(4),edge(i,j+1,k)(1),edge(i,j+1,k)(3),edge(i,j+1,k)(5),Phi(tet_grid.X(i,j,k))); // top left
            Append_Triangles(triangulated_surface,edge(i,j+1,k)(5),edge(i,j+1,k+1)(1),edge(i,j,k+1)(6),edge(i+1,j+1,k)(3),edge(i+1,j,k+1)(2),edge(i+1,j,k)(4),Phi(tet_grid.X(i,j+1,k+1))); // top right
            Append_Triangles(triangulated_surface,edge(i,j,k)(6),edge(i,j,k)(4),edge(i,j,k)(5),edge(i,j+1,k)(5),edge(i,j,k+1)(6),edge(i+1,j,k)(4),Phi(tet_grid.X(i,j,k)));} // center
    triangulated_surface.mesh.number_nodes=triangulated_surface.particles.array_collection->Size();
    triangulated_surface.Remove_Degenerate_Triangles();
}
//#####################################################################
// Function If_Zero_Crossing_Add_Particle
//#####################################################################
template<class T_GRID> int LEVELSET_3D<T_GRID>::
If_Zero_Crossing_Add_Particle(TRIANGULATED_SURFACE<T>& triangulated_surface,const TV& x1,const TV& x2) const
{
    int index=0;T phi1=Phi(x1),phi2=Phi(x2);
    if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
        index=triangulated_surface.particles.array_collection->Add_Element();
        triangulated_surface.particles.X(index)=LINEAR_INTERPOLATION<T,TV>::Linear(x1,x2,LEVELSET_UTILITIES<T>::Theta(phi1,phi2));}
    return index;
}
//#####################################################################
// Function Append_Triangles
//#####################################################################
// looking down at the node with phi1, e1-e2-e3 is conunter clockwise
// edges 1,2,3 come out of phi1 - edge4 is opposite edge3 - edge2 is opposite edge6 - edge1 is opposite edge5
template<class T_GRID> void LEVELSET_3D<T_GRID>::
Append_Triangles(TRIANGULATED_SURFACE<T>& triangulated_surface,const int e1,const int e2,const int e3,const int e4,const int e5,const int e6,const T phi1) const
{
    int number_positive=(e1>0)+(e2>0)+(e3>0)+(e4>0)+(e5>0)+(e6>0);if(number_positive == 0) return;assert(number_positive == 3 || number_positive == 4);
    if(e1 && e2 && e5 && e6){ // 2 triangles
        if(phi1 > 0){triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e1,e6,e2));triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e2,e6,e5));}
        else{triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e1,e2,e6));triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e2,e5,e6));}}
    else if(e2 && e3 && e4 && e6){
        if(phi1 > 0){triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e2,e4,e3));triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e3,e4,e6));}
        else{triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e2,e3,e4));triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e4,e3,e6));}}
    else if(e1 && e3 && e4 && e5){
        if(phi1 > 0){triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e1,e3,e5));triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e1,e5,e4));}
        else {triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e1,e5,e3));triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e1,e4,e5));}}
    else if(e1 && e2 && e3){ // 1 triangle
        if(phi1 > 0) triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e1,e3,e2));
        else triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e1,e2,e3));}
    else if(e1 && e4 && e6){
        if(phi1 > 0) triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e1,e6,e4));
        else triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e1,e4,e6));}
    else if(e3 && e5 && e6){
        if(phi1 > 0) triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e3,e5,e6));
        else triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e3,e6,e5));}
    else if (e4 && e2 && e5){
        if(phi1 > 0) triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e4,e5,e2));
        else triangulated_surface.mesh.elements.Append(VECTOR<int,3>(e4,e2,e5));}
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T_GRID> VECTOR<typename T_GRID::SCALAR,3> LEVELSET_3D<T_GRID>::
Normal(const VECTOR<T,3>& location) const
{
    if(normals) return normal_interpolation->Clamped_To_Array(grid,*normals,location).Normalized();
    else if(custom_normal_computation) return custom_normal_computation->Compute_Normal(*this,location);
    else return VECTOR<T,3>((Phi(VECTOR<T,3>(location.x+grid.dX.x,location.y,location.z))-Phi(VECTOR<T,3>(location.x-grid.dX.x,location.y,location.z)))/(2*grid.dX.x),
        (Phi(VECTOR<T,3>(location.x,location.y+grid.dX.y,location.z))-Phi(VECTOR<T,3>(location.x,location.y-grid.dX.y,location.z)))/(2*grid.dX.y),
        (Phi(VECTOR<T,3>(location.x,location.y,location.z+grid.dX.z))-Phi(VECTOR<T,3>(location.x,location.y,location.z-grid.dX.z)))/(2*grid.dX.z)).Normalized();
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class T_GRID> VECTOR<typename T_GRID::SCALAR,3> LEVELSET_3D<T_GRID>::
Extended_Normal(const VECTOR<T,3>& location) const
{
    if(normals) return normal_interpolation->Clamped_To_Array(grid,*normals,grid.Clamp(location)).Normalized();
    else return VECTOR<T,3>((Extended_Phi(VECTOR<T,3>(location.x+grid.dX.x,location.y,location.z))-Extended_Phi(VECTOR<T,3>(location.x-grid.dX.x,location.y,location.z)))/(2*grid.dX.x),
        (Extended_Phi(VECTOR<T,3>(location.x,location.y+grid.dX.y,location.z))-Extended_Phi(VECTOR<T,3>(location.x,location.y-grid.dX.y,location.z)))/(2*grid.dX.y),
        (Extended_Phi(VECTOR<T,3>(location.x,location.y,location.z+grid.dX.z))-Extended_Phi(VECTOR<T,3>(location.x,location.y,location.z-grid.dX.z)))/(2*grid.dX.z)).Normalized();
}
//#####################################################################
template class LEVELSET_3D<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_3D<GRID<VECTOR<double,3> > >;
#endif
