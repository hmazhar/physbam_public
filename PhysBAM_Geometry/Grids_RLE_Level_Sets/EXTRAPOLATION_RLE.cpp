//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/LINEAR_INTERPOLATION_RLE_HELPER.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/EXTRAPOLATION_RLE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> EXTRAPOLATION_RLE<T_GRID,T2>::
EXTRAPOLATION_RLE(const T_GRID& grid_input,ARRAY<T>& phi_input)
    :levelset(const_cast<T_GRID&>(grid_input),phi_input),use_default_value(false),neighbors_visible(0)
{
    Set_Band_Width();
    Set_Isobaric_Fix_Width();
    tolerance=levelset.small_number*levelset.grid.uniform_grid.dX.Min();
    VECTOR<T,3> DX(levelset.grid.uniform_grid.dX);
    dz_squared_over_dy_squared=sqr(DX.z/DX.y);dz_squared_over_dx_squared=sqr(DX.z/DX.x);dy_squared_over_dx_squared=sqr(DX.y/DX.x);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2> EXTRAPOLATION_RLE<T_GRID,T2>::
~EXTRAPOLATION_RLE()
{}
//#####################################################################
// Function Extrapolate_Cells
//#####################################################################
template<class T_GRID,class T2> void EXTRAPOLATION_RLE<T_GRID,T2>::
Extrapolate_Cells(ARRAY<T2>& u,const T time)
{
    ARRAY<T> phi_ghost(levelset.grid.number_of_cells,false);levelset.boundary->Fill_Ghost_Cells_Cell(levelset.grid,levelset.phi,phi_ghost,time);
    ARRAY<T2> u_ghost(levelset.grid.number_of_cells,false);boundary->Fill_Ghost_Cells_Cell(levelset.grid,u,u_ghost,time);
    ARRAY<bool,VECTOR<int,1> > done(0,levelset.grid.number_of_cells,false);
    done(0)=false;for(int n=1;n<=levelset.grid.number_of_cells;n++)done(n)=phi_ghost(n)<=-isobaric_fix_width;
    const ARRAY<VECTOR<int,T_GRID::number_of_neighbors_per_cell> >& neighbors=levelset.grid.Short_Cell_Neighbors();
    Extrapolate(neighbors,u_ghost,phi_ghost,done);
    ARRAY<T2>::Copy(u_ghost,u);boundary->Apply_Boundary_Condition(levelset.grid,u,time);
}
//#####################################################################
// Function Extrapolate_Faces
//#####################################################################
template<class T_GRID,class T2> void EXTRAPOLATION_RLE<T_GRID,T2>::
Extrapolate_Faces(ARRAY<T2>& u,const ARRAY<bool>& fixed) const
{
    ARRAY<T> phi_ghost_face(levelset.grid.number_of_faces);LINEAR_INTERPOLATION_RLE_HELPER<T_GRID>::Interpolate_From_Short_Cells_To_Short_Faces(levelset.grid,levelset.phi,phi_ghost_face);
    ARRAY<bool,VECTOR<int,1> > done(0,levelset.grid.number_of_faces,false);done(0)=false;
    for(int n=1;n<=levelset.grid.number_of_faces;n++)done(n)=fixed(n);
    const ARRAY<VECTOR<int,T_GRID::number_of_neighbors_per_cell> >& neighbors=levelset.grid.Short_Face_Neighbors();
    Extrapolate(neighbors,u,phi_ghost_face,done);
}
//#####################################################################
// Function Extrapolate
//#####################################################################
template<class T_GRID,class T2> void EXTRAPOLATION_RLE<T_GRID,T2>::
Extrapolate(const T_FACE_NEIGHBORS& neighbors,ARRAY<T2>& u,const ARRAY<T>& phi,ARRAY<bool,VECTOR<int,1> >& done) const
{
    assert(!neighbors_visible || neighbors.m==neighbors_visible->m);
    if(use_default_value) for(int n=1;n<=phi.m;n++)if(!done(n)) u(n)=default_value;
    int heap_length=0;
    ARRAY<bool> close(phi.m);
    ARRAY<int> heap(phi.m);
    for(int n=1;n<=phi.m;n++)if(done(n))for(int i=1;i<=T_GRID::number_of_neighbors_per_cell;i++){int neighbor=neighbors(n)(i);
        if(Neighbor_Visible(i,n,neighbor) && !done(neighbor) && !close(neighbor))
            Add_To_Heap(phi,heap,heap_length,close,neighbor);}

    while(heap_length && phi(heap(1))<=band_width){
        int node=Remove_Root_From_Heap(phi,heap,heap_length,close);
        done(node)=true;
        Update_Close_Point(neighbors,u,phi,done,node);
        for(int i=1;i<=T_GRID::number_of_neighbors_per_cell;i++){int neighbor=neighbors(node)(i);
            if(Neighbor_Visible(i,node,neighbor) && !done(neighbor) && !close(neighbor))
                Add_To_Heap(phi,heap,heap_length,close,neighbor);}}
}
//#####################################################################
// Function Update_Close_Point
//#####################################################################
template<class T_GRID,class T2> inline void EXTRAPOLATION_RLE<T_GRID,T2>::
Update_Close_Point(const ARRAY<VECTOR<int,4> >& neighbors,ARRAY<T2>& u,const ARRAY<T>& phi,const ARRAY<bool,VECTOR<int,1> >& done,const int node) const
{
    int n_mx,n_px,n_my,n_py;neighbors(node).Get(n_mx,n_px,n_my,n_py);
    bool check_left=done(n_mx),check_right=done(n_px),check_bottom=done(n_my),check_top=done(n_py);
    if(neighbors_visible){
        if(check_left && !(*neighbors_visible)(n_mx)(1)) check_left=false;
        if(check_right && !(*neighbors_visible)(node)(1)) check_right=false;
        if(check_bottom && !(*neighbors_visible)(n_my)(2)) check_bottom=false;
        if(check_top && !(*neighbors_visible)(node)(2)) check_top=false;}
    T2 value_x=T2(),value_y=T2();T phix_dx=0,phiy_dy=0;
    bool no_x=!check_left && !check_right;if(!no_x){int nx=check_left?(check_right?(phi(n_mx)<=phi(n_px)?n_mx:n_px):n_mx):n_px;value_x=u(nx);phix_dx=phi(node)-phi(nx);}
    bool no_y=!check_bottom && !check_top;if(!no_y){int ny=check_bottom?(check_top?(phi(n_my)<=phi(n_py)?n_my:n_py):n_my):n_py;value_y=u(ny);phiy_dy=phi(node)-phi(ny);}
    u(node)=Solve_Close_Point(no_x,value_x,phix_dx,no_y,value_y,phiy_dy,dy_squared_over_dx_squared,tolerance);
}
//#####################################################################
// Function Update_Close_Point
//#####################################################################
template<class T_GRID,class T2> inline void EXTRAPOLATION_RLE<T_GRID,T2>::
Update_Close_Point(const ARRAY<VECTOR<int,6> >& neighbors,ARRAY<T2>& u,const ARRAY<T>& phi,const ARRAY<bool,VECTOR<int,1> >& done,const int node) const
{
    int n_mx,n_px,n_my,n_py,n_mz,n_pz;neighbors(node).Get(n_mx,n_px,n_my,n_py,n_mz,n_pz);
    bool check_left=done(n_mx),check_right=done(n_px),check_bottom=done(n_my),check_top=done(n_py),check_front=done(n_mz),check_back=done(n_pz);
    if(neighbors_visible){
        if(check_left && !(*neighbors_visible)(n_mx)(1)) check_left=false;
        if(check_right && !(*neighbors_visible)(node)(1)) check_right=false;
        if(check_bottom && !(*neighbors_visible)(n_my)(2)) check_bottom=false;
        if(check_top && !(*neighbors_visible)(node)(2)) check_top=false;
        if(check_front && !(*neighbors_visible)(n_mz)(3)) check_front=false;
        if(check_back && !(*neighbors_visible)(node)(3)) check_back=false;}
    T2 value_x=T2(),value_y=T2(),value_z=T2();T phix_dx=0,phiy_dy=0,phiz_dz=0;
    bool no_x=!check_left && !check_right;if(!no_x){int nx=check_left?(check_right?(phi(n_mx)<=phi(n_px)?n_mx:n_px):n_mx):n_px;value_x=u(nx);phix_dx=phi(node)-phi(nx);}
    bool no_y=!check_bottom && !check_top;if(!no_y){int ny=check_bottom?(check_top?(phi(n_my)<=phi(n_py)?n_my:n_py):n_my):n_py;value_y=u(ny);phiy_dy=phi(node)-phi(ny);}
    bool no_z=!check_front && !check_back;if(!no_z){int nz=check_front?(check_back?(phi(n_mz)<=phi(n_pz)?n_mz:n_pz):n_mz):n_pz;value_z=u(nz);phiz_dz=phi(node)-phi(nz);}
    u(node)=Solve_Close_Point(no_x,value_x,phix_dx,no_y,value_y,phiy_dy,no_z,value_z,phiz_dz,dz_squared_over_dy_squared,dz_squared_over_dx_squared,dy_squared_over_dx_squared,tolerance);
}
//#####################################################################
template class EXTRAPOLATION_RLE<RLE_GRID_2D<float>,float>;
template class EXTRAPOLATION_RLE<RLE_GRID_3D<float>,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EXTRAPOLATION_RLE<RLE_GRID_2D<double>,double>;
template class EXTRAPOLATION_RLE<RLE_GRID_3D<double>,double>;
#endif
#endif
