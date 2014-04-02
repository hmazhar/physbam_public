#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> BOUNDARY_MPI<T_GRID,T2>::
BOUNDARY_MPI(T_MPI_GRID* mpi_grid_input,T_BOUNDARY_T2& boundary_input)
    :mpi_grid(mpi_grid_input),boundary(boundary_input)
{
    assert(&boundary);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2> BOUNDARY_MPI<T_GRID,T2>::
~BOUNDARY_MPI()
{
}
//#####################################################################
// Function Set_Constant_Extrapolation
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_MPI<T_GRID,T2>::
Set_Constant_Extrapolation(const TV_SIDES& constant_extrapolation_input)
{
    boundary.Set_Constant_Extrapolation(constant_extrapolation_input);
}
//#####################################################################
// Function Constant_Extrapolation
//#####################################################################
template<class T_GRID,class T2> bool BOUNDARY_MPI<T_GRID,T2>::
Constant_Extrapolation(const int side) const
{
    return boundary.Constant_Extrapolation(side);
}
//#####################################################################
// Function Set_Fixed_Boundary
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_MPI<T_GRID,T2>::
Set_Fixed_Boundary(const bool use_fixed_boundary_input,const T2 fixed_boundary_value_input)
{
    boundary.Set_Fixed_Boundary(use_fixed_boundary_input,fixed_boundary_value_input);
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_MPI<T_GRID,T2>::
Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells_input)
{
    boundary.Fill_Ghost_Cells(grid,u,u_ghost,dt,time,number_of_ghost_cells_input);
    mpi_grid->Exchange_Boundary_Cell_Data(u_ghost,number_of_ghost_cells_input);
}
//#####################################################################
// Function Fill_Ghost_Cells_Face
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_MPI<T_GRID,T2>::
Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells_input)
{
    boundary.Fill_Ghost_Cells_Face(grid,u,u_ghost,time,number_of_ghost_cells_input);
    mpi_grid->Exchange_Boundary_Face_Data(u_ghost,number_of_ghost_cells_input);
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_MPI<T_GRID,T2>::
Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_T2& u,const T time)
{
    boundary.Apply_Boundary_Condition(grid,u,time);
}
//#####################################################################
// Function Apply_Boundary_Condition_Face
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_MPI<T_GRID,T2>::
Apply_Boundary_Condition_Face(const T_GRID& grid,T_FACE_ARRAYS_T2& u,const T time)
{
    boundary.Apply_Boundary_Condition_Face(grid,u,time);
    mpi_grid->Average_Common_Face_Data(u);
}
template class BOUNDARY_MPI<GRID<VECTOR<float,1> >,VECTOR<float,1> >;
template class BOUNDARY_MPI<GRID<VECTOR<float,1> >,VECTOR<float,3> >;
template class BOUNDARY_MPI<GRID<VECTOR<float,1> >,float>;
template class BOUNDARY_MPI<GRID<VECTOR<float,2> >,VECTOR<float,2> >;
template class BOUNDARY_MPI<GRID<VECTOR<float,2> >,VECTOR<float,4> >;
template class BOUNDARY_MPI<GRID<VECTOR<float,2> >,float>;
template class BOUNDARY_MPI<GRID<VECTOR<float,3> >,VECTOR<float,3> >;
template class BOUNDARY_MPI<GRID<VECTOR<float,3> >,VECTOR<float,5> >;
template class BOUNDARY_MPI<GRID<VECTOR<float,3> >,float>;
template class BOUNDARY_MPI<GRID<VECTOR<float,1> >,MATRIX<float,1,1> >;
template class BOUNDARY_MPI<GRID<VECTOR<float,2> >,SYMMETRIC_MATRIX<float,2> >;
template class BOUNDARY_MPI<GRID<VECTOR<float,3> >,SYMMETRIC_MATRIX<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BOUNDARY_MPI<GRID<VECTOR<double,1> >,VECTOR<double,1> >;
template class BOUNDARY_MPI<GRID<VECTOR<double,1> >,VECTOR<double,3> >;
template class BOUNDARY_MPI<GRID<VECTOR<double,1> >,double>;
template class BOUNDARY_MPI<GRID<VECTOR<double,2> >,VECTOR<double,2> >;
template class BOUNDARY_MPI<GRID<VECTOR<double,2> >,VECTOR<double,4> >;
template class BOUNDARY_MPI<GRID<VECTOR<double,2> >,double>;
template class BOUNDARY_MPI<GRID<VECTOR<double,3> >,VECTOR<double,3> >;
template class BOUNDARY_MPI<GRID<VECTOR<double,3> >,VECTOR<double,5> >;
template class BOUNDARY_MPI<GRID<VECTOR<double,3> >,double>;
template class BOUNDARY_MPI<GRID<VECTOR<double,1> >,MATRIX<double,1,1> >;
template class BOUNDARY_MPI<GRID<VECTOR<double,2> >,SYMMETRIC_MATRIX<double,2> >;
template class BOUNDARY_MPI<GRID<VECTOR<double,3> >,SYMMETRIC_MATRIX<double,3> >;
#endif
