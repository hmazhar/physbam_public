//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CATMULL_ROM_SPLINE_INTERPOLATION.h
//#####################################################################
#ifndef __CATMULL_ROM_SPLINE_INTERPOLATION__
#define __CATMULL_ROM_SPLINE_INTERPOLATION__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<T_GRID>
class CATMULL_ROM_SPLINE_INTERPOLATION:public INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T::template REBIND<bool>::TYPE TV_BOOL;
    typedef typename T_GRID::ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;typedef typename T_GRID::CELL_ITERATOR T_CELL_ITERATOR;
public:
    typedef INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP> BASE;
    using BASE::Clamped_Index_End_Minus_One;

    T tension;

    CATMULL_ROM_SPLINE_INTERPOLATION()
        :tension((T).5)
    {}

    T2 Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE
    {return From_Base_Node(grid,u,X,Clamped_Index_End_Minus_One(grid,u,X));}

    T2 Clamped_To_Array_Derivative(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_BOOL& derivatives) const
    {return From_Base_Node_Derivative(grid,u,X,Clamped_Index_End_Minus_One(grid,u,X),derivatives);}

    T2 From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& index) const PHYSBAM_OVERRIDE
    {T basis[T_GRID::dimension][4];T2 sum=T2();TV X_normalized=(X-grid.X(index))*grid.one_over_dX;
    for(int axis=1;axis<=T_GRID::dimension;axis++) Catmull_Rom_Basis(X_normalized[axis],basis[axis-1]);
    for(T_CELL_ITERATOR iterator(grid,RANGE<TV_INT>(index-TV_INT::All_Ones_Vector(),index+2*TV_INT::All_Ones_Vector()));iterator.Valid();iterator.Next()){
        T product=1;for(int axis=1;axis<=T_GRID::dimension;axis++) product*=basis[axis-1][iterator.Cell_Index()[axis]-index[axis]+1];
        sum+=u(iterator.Cell_Index())*product;}
    return sum;}

    T2 From_Base_Node_Derivative(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& index,const TV_BOOL& derivative) const
    {T basis[T_GRID::dimension][4];T2 sum=T2();TV X_normalized=(X-grid.X(index))*grid.one_over_dX;
    for(int axis=1;axis<=T_GRID::dimension;axis++)
        if(derivative[axis]) Catmull_Rom_Basis_Derivative(X_normalized[axis],basis[axis-1]);
        else Catmull_Rom_Basis(X_normalized[axis],basis[axis-1]);
    for(T_CELL_ITERATOR iterator(grid,RANGE<TV_INT>(index-TV_INT::All_Ones_Vector(),index+2*TV_INT::All_Ones_Vector()));iterator.Valid();iterator.Next()){
        T product=1;for(int axis=1;axis<=T_GRID::dimension;axis++) product*=basis[axis-1][iterator.Cell_Index()[axis]-index[axis]+1];
        sum+=u(iterator.Cell_Index())*product;}
    return sum;}

private:
    void Catmull_Rom_Basis(const T u,T* basis) const
    {T sqr_u=sqr(u),cube_u=cube(u);
    basis[0]=tension*(-u+2*sqr_u-cube_u);basis[1]=1+(tension-3)*sqr_u+(2-tension)*cube_u;basis[2]=tension*u+(3-2*tension)*sqr_u+(tension-2)*cube_u;basis[3]=tension*(-sqr_u+cube_u);}

    void Catmull_Rom_Basis_Derivative(const T u,T* basis) const
    {T sqr_u=sqr(u);
    basis[0]=tension*(-1+4*u-3*sqr_u);basis[1]=(2*tension-6)*u+(6-3*tension)*sqr_u;basis[2]=tension+(6-4*tension)*u+(3*tension-6)*sqr_u;basis[3]=tension*(-2*u+3*sqr_u);}

//#####################################################################
};
}
#endif
