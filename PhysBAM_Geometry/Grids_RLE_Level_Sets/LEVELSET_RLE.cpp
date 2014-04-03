//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_TRANSFER_ITERATOR.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/FAST_MARCHING_METHOD_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
template<class T,class T_GRID> typename LEVELSET<T,T_GRID>::T_LINEAR_INTERPOLATION_SCALAR LEVELSET<T,T_GRID>::interpolation_default;
template<class T,class T_GRID> typename REBIND<typename LEVELSET<T,T_GRID>::T_LINEAR_INTERPOLATION_SCALAR,typename T_GRID::VECTOR_T>::TYPE LEVELSET<T,T_GRID>::normal_interpolation_default;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LEVELSET_RLE<T_GRID>::
LEVELSET_RLE(T_GRID& grid_input,ARRAY<T>& phi_input)
    :grid(grid_input),phi(phi_input),normals(0),curvature(0),block_minimum(0),block_maximum(0)
{   
    Set_Band_Width(20);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> LEVELSET_RLE<T_GRID>::
~LEVELSET_RLE()
{
    Clean_Memory();
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T_GRID> void LEVELSET_RLE<T_GRID>::
Clean_Memory()
{
    delete normals;normals=0;delete curvature;curvature=0;delete block_minimum;block_minimum=0;delete block_maximum;block_maximum=0;
}
//#####################################################################
// Function Collision_Aware_Phi
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_RLE<T_GRID>::
Collision_Aware_Phi(const T_BLOCK& block,const TV& location) const
{
    assert(collision_body_list);
    return LINEAR_INTERPOLATION_COLLIDABLE_CELL_RLE<T_GRID,T>(*collision_body_list,&valid_mask_current,collidable_phi_replacement_value).From_Block_Cell(block,phi,location);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T LEVELSET_RLE<T_GRID>::
Normal(const TV& X) const
{   
    if(normals){
        T_BLOCK block(grid,X);if(!block) return TV::Axis_Vector(1);
        return Normal(block,X);}
    TV normal;
    for(int axis=1;axis<=T_GRID::dimension;axis++){TV offset=grid.uniform_grid.dX[axis]*TV::Axis_Vector(axis);
        normal[axis]=(Phi(X+offset)-Phi(X-offset))*grid.uniform_grid.one_over_dX[axis];}
    return normal.Normalized();
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T LEVELSET_RLE<T_GRID>::
Extended_Normal(const TV& X) const
{   
    if(normals) return LINEAR_INTERPOLATION_RLE_HELPER<T_GRID>::Interpolate_Cells_Short(grid,*normals,grid.Clamp(X));
    TV normal;
    for(int axis=1;axis<=T_GRID::dimension;axis++){TV offset=grid.uniform_grid.dX[axis]*TV::Axis_Vector(axis);
        normal[axis]=(Extended_Phi(X+offset)-Extended_Phi(X-offset))*grid.uniform_grid.one_over_dX[axis];}
    return normal.Normalized();
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_RLE<T_GRID>::
Phi_Secondary(const TV& X) const
{
    T_BLOCK block(grid,X);if(block) return secondary_interpolation->From_Block_Cell(block,phi,X);
    CELL_ITERATOR cell(grid,X);return phi(cell.Cell());
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_RLE<T_GRID>::
CFL(const ARRAY<T>& V) const
{
    T max_V_norm=0;
    for(CELL_ITERATOR cell(grid,0);cell;cell++){
        T local_V_norm=0;
        for(int axis=1;axis<=T_GRID::dimension;axis++)local_V_norm+=maxabs(V(cell.First_Face_Index(axis)),V(cell.Second_Face_Index(axis)));
        max_V_norm=max(max_V_norm,local_V_norm);}
    T dt_convection=max_V_norm/grid.uniform_grid.dX.x;
    T dt_curvature=curvature_motion?sigma*T_GRID::number_of_faces_per_cell/sqr(grid.uniform_grid.dX.x):0;
    T dt_overall=dt_convection+dt_curvature;
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function Curvature_CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_RLE<T_GRID>::
Curvature_CFL() const
{
    T dt_curvature=curvature_motion?sigma*T_GRID::number_of_faces_per_cell/sqr(grid.uniform_grid.dX.x):0;
    return Robust_Divide((T)1,dt_curvature);
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class T_GRID> void LEVELSET_RLE<T_GRID>::
Compute_Normals(const T time)
{
    ARRAY<T> phi_ghost(grid.number_of_cells,false);boundary->Fill_Ghost_Cells_Cell(grid,phi,phi_ghost,time);
    Compute_Normals(phi_ghost);
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class T_GRID> struct Compute_Normals_Helper{template<class T_FACE> static void Apply(const T_GRID& grid,const ARRAY<typename T_GRID::SCALAR>& phi,ARRAY<typename T_GRID::VECTOR_T>& normals)
{
    int axis=T_FACE::Axis();
    T_FACE face(grid,3,false);
    if(!face.cell1.run->cell){LOG::cout<<"Nasty bug detected.  There's some chance this is a gcc optimization bug."<<std::endl;PHYSBAM_FATAL_ERROR();}
    for(;face;face++)if(face.Both_Cells_Short()){
        int c1=face.cell1.Cell(),c2=face.cell2.Cell();normals(c2)[axis]=-phi(c1);normals(c1)[axis]+=phi(c2);}
}};
template<class T_GRID> void LEVELSET_RLE<T_GRID>::
Compute_Normals(ARRAY<T>& phi_ghost)
{
    if(!normals) normals=new ARRAY<TV>(grid.number_of_cells);else normals->Resize(grid.number_of_cells,false,false);
    T_GRID::template Face_Loop<Compute_Normals_Helper<T_GRID> >(grid,phi_ghost,*normals);
    TV dx_over_DX=grid.uniform_grid.dX.x*grid.uniform_grid.one_over_dX;
    for(int c=1;c<=grid.number_of_cells;c++){
        for(int axis=2;axis<=T_GRID::dimension;axis++)(*normals)(c)[axis]*=dx_over_DX[axis];
        (*normals)(c).Normalize();}
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T_GRID> void LEVELSET_RLE<T_GRID>::
Compute_Curvature(const T time)
{
    ARRAY<T> phi_ghost(grid.number_of_cells,false);boundary->Fill_Ghost_Cells_Cell(grid,phi,phi_ghost,time);
    Compute_Curvature(phi_ghost);
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class T_GRID> void LEVELSET_RLE<T_GRID>::
Compute_Curvature(const ARRAY<T>& phi_ghost)
{
    if(!curvature) curvature=new ARRAY<T>(grid.number_of_cells);else curvature->Resize(grid.number_of_cells);
    Compute_Curvature_Helper(grid,phi_ghost,small_number,*curvature);
}
template<class T> static void Compute_Curvature_Helper(const RLE_GRID_2D<T>& grid,const ARRAY<T>& phi,const T small_number,ARRAY<T>& curvature)
{
    T dx=grid.uniform_grid.dX.x,dy=grid.uniform_grid.dX.y,one_over_two_dx=1/(2*dx),one_over_two_dy=1/(2*dy),one_over_dx_squared=1/sqr(dx),one_over_dy_squared=1/sqr(dy),
        one_over_four_dx_dy=1/(4*dx*dy),max_curvature=1/grid.uniform_grid.dX.Max();
    for(typename RLE_GRID_2D<T>::CELL_ITERATOR cell(grid,2);cell;cell++){
        int c=cell.Cell(),i=cell.i,j=cell.j;
        T phi_c=phi(c),phi_mx=grid.Cell_Value(phi,i-1,j),phi_px=grid.Cell_Value(phi,i+1,j),phi_my=grid.Cell_Value(phi,i,j-1),
            phi_py=grid.Cell_Value(phi,i,j+1),phi_mx_my=grid.Cell_Value(phi,i-1,j-1),phi_px_my=grid.Cell_Value(phi,i+1,j-1),
            phi_mx_py=grid.Cell_Value(phi,i-1,j+1),phi_px_py=grid.Cell_Value(phi,i+1,j+1); // TODO: optimize
        T phix=(phi_px-phi_mx)*one_over_two_dx,phiy=(phi_py-phi_my)*one_over_two_dy,phixx=(phi_px-2*phi_c+phi_mx)*one_over_dx_squared,
            phiyy=(phi_py-2*phi_c+phi_my)*one_over_dy_squared,phixy=(phi_px_py-phi_px_my-phi_mx_py+phi_mx_my)*one_over_four_dx_dy;
        T denominator=sqrt(sqr(phix)+sqr(phiy)),curvature_value;
        if(denominator >= small_number) curvature_value=-(sqr(phix)*phiyy-2*phix*phiy*phixy+sqr(phiy)*phixx)/cube(denominator);
        else curvature_value=LEVELSET_UTILITIES<T>::Sign(phi_c)*max_curvature;
        curvature(c)=minmag(curvature_value,LEVELSET_UTILITIES<T>::Sign(curvature_value)*max_curvature);}
}
template<class T> static void Compute_Curvature_Helper(const RLE_GRID_3D<T>& grid,const ARRAY<T>& phi,const T small_number,ARRAY<T>& curvature)
{
    T dx=grid.uniform_grid.dX.x,dy=grid.uniform_grid.dX.y,dz=grid.uniform_grid.dX.z,one_over_two_dx=1/(2*dx),one_over_two_dy=1/(2*dy),one_over_two_dz=1/(2*dz),one_over_dx_squared=1/sqr(dx),
        one_over_dy_squared=1/sqr(dy),one_over_dz_squared=1/sqr(dz),one_over_four_dx_dy=1/(4*dx*dy),one_over_four_dx_dz=1/(4*dx*dz),one_over_four_dy_dz=1/(4*dy*dz),
        max_curvature=1/grid.uniform_grid.dX.Max();
    for(typename RLE_GRID_3D<T>::CELL_ITERATOR cell(grid,2);cell;cell++){
        int c=cell.Cell(),i=cell.i,j=cell.j,ij=cell.ij;
        T phi_c=phi(c),phi_mx=grid.Cell_Value(phi,i-1,j,ij),phi_px=grid.Cell_Value(phi,i+1,j,ij),phi_my=grid.Cell_Value(phi,i,j-1,ij),
            phi_py=grid.Cell_Value(phi,i,j+1,ij),phi_mz=grid.Cell_Value(phi,i,j,ij-1),phi_pz=grid.Cell_Value(phi,i,j,ij+1),
            phi_mx_my=grid.Cell_Value(phi,i-1,j-1,ij),phi_px_my=grid.Cell_Value(phi,i+1,j-1,ij),
            phi_mx_py=grid.Cell_Value(phi,i-1,j+1,ij),phi_px_py=grid.Cell_Value(phi,i+1,j+1,ij),
            phi_mx_mz=grid.Cell_Value(phi,i-1,j,ij-1),phi_px_mz=grid.Cell_Value(phi,i+1,j,ij-1),
            phi_mx_pz=grid.Cell_Value(phi,i-1,j,ij+1),phi_px_pz=grid.Cell_Value(phi,i+1,j,ij+1),
            phi_my_mz=grid.Cell_Value(phi,i,j-1,ij-1),phi_py_mz=grid.Cell_Value(phi,i,j+1,ij-1),
            phi_my_pz=grid.Cell_Value(phi,i,j-1,ij+1),phi_py_pz=grid.Cell_Value(phi,i,j+1,ij+1); // TODO: optimize
        T phix=(phi_px-phi_mx)*one_over_two_dx,phiy=(phi_py-phi_my)*one_over_two_dy,phiz=(phi_pz-phi_mz)*one_over_two_dz,
            phixx=(phi_px-2*phi_c+phi_mx)*one_over_dx_squared,phiyy=(phi_py-2*phi_c+phi_my)*one_over_dy_squared,phizz=(phi_pz-2*phi_c+phi_mz)*one_over_dz_squared,
            phixy=(phi_px_py-phi_px_my-phi_mx_py+phi_mx_my)*one_over_four_dx_dy,phixz=(phi_px_pz-phi_px_mz-phi_mx_pz+phi_mx_mz)*one_over_four_dx_dz,
            phiyz=(phi_py_pz-phi_py_mz-phi_my_pz+phi_my_mz)*one_over_four_dy_dz;
        T denominator=sqrt(sqr(phix)+sqr(phiy)+sqr(phiz));
        if(denominator >= small_number){
            T curvature_value=-(sqr(phix)*phiyy-2*phix*phiy*phixy+sqr(phiy)*phixx+sqr(phix)*phizz-2*phix*phiz*phixz+sqr(phiz)*phixx
                +sqr(phiy)*phizz-2*phiy*phiz*phiyz+sqr(phiz)*phiyy)/cube(denominator);
            curvature(c)=minmag(curvature_value,LEVELSET_UTILITIES<T>::Sign(curvature_value)*max_curvature);}
        else curvature(c)=LEVELSET_UTILITIES<T>::Sign(phi_c)*max_curvature;}
}
//#####################################################################
// Function Compute_Block_Minimum_And_Maximum
//#####################################################################
template<class T_GRID> void LEVELSET_RLE<T_GRID>::
Compute_Block_Minimum_And_Maximum(const bool recompute_if_exists)
{
    if(!recompute_if_exists && block_minimum && block_minimum) return;
    if(!block_minimum) block_minimum=new ARRAY<T>;block_minimum->Resize(grid.number_of_blocks,false,false);
    if(!block_maximum) block_maximum=new ARRAY<T>;block_maximum->Resize(grid.number_of_blocks,false,false);
    for(BLOCK_ITERATOR block(grid,grid.number_of_ghost_cells);block;block++){int b=block.Block();
        int c1=block.Cell(0),c2=block.Cell(1);
        T min_phi,max_phi;
        if(phi(c1)<phi(c2)){min_phi=phi(c1);max_phi=phi(c2);}else{min_phi=phi(c2);max_phi=phi(c1);}
        for(int i=2;i<T_GRID::number_of_cells_per_block;i++){
            T phi_c=phi(block.Cell(i));if(phi_c<min_phi) min_phi=phi_c;else if(phi_c>max_phi) max_phi=phi_c;}
        (*block_minimum)(b)=min_phi;(*block_maximum)(b)=max_phi;}
}
//#####################################################################
// Function Curvature_Motion
//#####################################################################
template<class T_GRID> void LEVELSET_RLE<T_GRID>::
Curvature_Motion(const T dt,const T time)
{
    ARRAY<T> phi_ghost(grid.number_of_cells,false);boundary->Fill_Ghost_Cells_Cell(grid,phi,phi_ghost,time);
    Curvature_Motion(dt,phi_ghost);
}
//#####################################################################
// Function Curvature_Motion
//#####################################################################
template<class T_GRID> void LEVELSET_RLE<T_GRID>::
Curvature_Motion(const T dt,const ARRAY<T>& phi_ghost)
{
    TV one_over_two_DX=(T).5*grid.uniform_grid.one_over_dX;
    bool curvature_defined=curvature!=0;Compute_Curvature(phi_ghost);
    for(CELL_ITERATOR cell(grid,0);cell;cell++)if(cell.Short()){int c=cell.Cell();TV_INT I=cell.I();
        T gradient_magnitude_squared=0;
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT offset=TV_INT::Axis_Vector(axis);
            gradient_magnitude_squared+=sqr(one_over_two_DX[axis]*(grid.Cell_Value(phi_ghost,I+offset)-grid.Cell_Value(phi_ghost,I-offset)));}
        phi(c)-=dt*sigma*(*curvature)(c)*sqrt(gradient_magnitude_squared);} 
    if(!curvature_defined){delete curvature;curvature=0;}
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T_GRID> void LEVELSET_RLE<T_GRID>::
Fast_Marching_Method(const T time)
{
    boundary->Fill_Ghost_Cells_Cell(grid,phi,phi,time);
    T stopping_distance=max(half_bandwidth,grid.negative_bandwidth,grid.positive_bandwidth)+grid.uniform_grid.dX.Max();
    FAST_MARCHING_METHOD_RLE<T_GRID> fast_marching(*this);
    Compute_Normals(phi);
    fast_marching.Slow_Marching(phi,stopping_distance); // TODO: somewhat odd for Fast_Marching_Method to unconditionally call Slow_Marching
    boundary->Apply_Boundary_Condition(grid,phi,time);
}
//#####################################################################
// Function Approximate_Negative_Size
//#####################################################################
// calculates the approximate volume using Heaviside functions
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_RLE<T_GRID>::
Approximate_Negative_Size(const T interface_thickness,const T time) const
{
    T interface_half_width=interface_thickness*grid.uniform_grid.dX.Max()/2,volume=0;
    for(CELL_ITERATOR cell(grid,0);cell;cell++)volume+=cell.length*LEVELSET_UTILITIES<T>::Heaviside(-phi(cell.Cell()),interface_half_width);
    return volume*grid.uniform_grid.Cell_Size();
}
//#####################################################################
// Function Approximate_Positive_Size
//#####################################################################
// calculates the approximate volume using Heaviside functions
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_RLE<T_GRID>::
Approximate_Positive_Size(const T interface_thickness,const T time) const
{
    T interface_half_width=interface_thickness*grid.uniform_grid.dX.Max()/2,volume=0;
    for(CELL_ITERATOR cell(grid,0);cell;cell++)volume+=cell.length*LEVELSET_UTILITIES<T>::Heaviside(phi(cell.Cell()),interface_half_width);
    return volume*grid.uniform_grid.Cell_Size();
}
//#####################################################################
// Function Transfer_Phi
//#####################################################################
template<class T_GRID> void LEVELSET_RLE<T_GRID>::
Transfer_Phi(const T_GRID& new_grid)
{
    ARRAY<T> new_phi(new_grid.number_of_cells);
    if(valid_mask_current.m){
        ARRAY<bool> new_valid_mask(new_grid.number_of_cells,false);ARRAYS_COMPUTATIONS::Fill(new_valid_mask,true);
        for(RLE_GRID_TRANSFER_ITERATOR<T_GRID,CELL_ITERATOR> cells(grid,new_grid,grid.number_of_ghost_cells);cells;cells++){
            Transfer_Cells(cells,phi,new_phi);
            if(cells.source.Short() && cells.destination.Short()) new_valid_mask(cells.destination.Cell())=valid_mask_current(cells.source.Cell());}
        ARRAY<bool>::Exchange_Arrays(valid_mask_current,new_valid_mask);}
    else for(RLE_GRID_TRANSFER_ITERATOR<T_GRID,CELL_ITERATOR> cells(grid,new_grid,grid.number_of_ghost_cells);cells;cells++)Transfer_Cells(cells,phi,new_phi);
    ARRAY<T>::Exchange_Arrays(phi,new_phi);
}
//#####################################################################
// Function Transfer_Phi
//#####################################################################
template<class T_GRID> void LEVELSET_RLE<T_GRID>::
Rebuild_Grid(const int extra_cells,const T_ARRAYS_HORIZONTAL_INT* ground_j)
{
    ARRAY<bool> cell_should_be_long(grid.number_of_cells);
    for(int c=1;c<=grid.number_of_cells;c++)cell_should_be_long(c)=phi(c)>=grid.positive_bandwidth || phi(c)<=-grid.negative_bandwidth;
    T_GRID new_grid;
    T_GRID::Rebuild(grid,new_grid,cell_should_be_long,extra_cells,ground_j);
    Transfer_Phi(new_grid);
    T_GRID::Transfer(new_grid,grid);
}
//#####################################################################
template class LEVELSET_RLE<RLE_GRID_2D<float> >;
template class LEVELSET_RLE<RLE_GRID_3D<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_RLE<RLE_GRID_2D<double> >;
template class LEVELSET_RLE<RLE_GRID_3D<double> >;
#endif
#endif
