//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_MAC_3D_HELPER
//#####################################################################
#ifndef __LINEAR_INTERPOLATION_MAC_3D_HELPER__
#define __LINEAR_INTERPOLATION_MAC_3D_HELPER__

#include <PhysBAM_Tools/Grids_Uniform/BLOCK_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

template<class T_GRID> class DYADIC_GRID_ITERATOR_FACE;

template<class T_GRID>
class LINEAR_INTERPOLATION_MAC_3D_HELPER
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::BLOCK T_BLOCK;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS;typedef typename T_FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS;typedef typename T_GRID::INDEX T_INDEX;
private:
    TV base,center,one_over_DX;
    T u2,u5,u8,u11,slope_u12,slope_u23,slope_u45,slope_u56,slope_u78,slope_u89,slope_u10_11,slope_u11_12; // standard z-y-x major ordering
    T v2,v5,v8,v11,slope_v12,slope_v23,slope_v45,slope_v56,slope_v78,slope_v89,slope_v10_11,slope_v11_12; // z-x-y major ordering for symmetry with u
    T w2,w5,w8,w11,slope_w12,slope_w23,slope_w45,slope_w56,slope_w78,slope_w89,slope_w10_11,slope_w11_12; // y-x-z major ordering for symmetry with u
public:

    LINEAR_INTERPOLATION_MAC_3D_HELPER(const T_BLOCK& block,const T_FACE_ARRAYS& face_velocities);
    ~LINEAR_INTERPOLATION_MAC_3D_HELPER();

    template<class T_FACE_LOOKUP>
    static TV Interpolate_Face(const T_BLOCK& block,const T_FACE_LOOKUP& face_velocities,const TV& X)
    {return Interpolate_Face_Transformed(block,face_velocities,Transformed(block,X));}

    template<class T_FACE_LOOKUP>
    static T Interpolate_Face_Component(const int axis,const T_BLOCK& block,const T_FACE_LOOKUP& face_velocities,const TV& X)
    {switch(axis){
      case 1:return Interpolate_Face_X_Transformed(block,face_velocities,Transformed(block,X));
      case 2:return Interpolate_Face_Y_Transformed(block,face_velocities,Transformed(block,X));
      default:assert(axis==3);return Interpolate_Face_Z_Transformed(block,face_velocities,Transformed(block,X));}}

    template<class T_FACE_LOOKUP>
    static ARRAY<PAIR<FACE_INDEX<3>,T> > Interpolate_Face_Component_Weights(const int axis,const T_BLOCK& block,const T_FACE_LOOKUP& face_velocities,const TV& X)
    {switch(axis){
      case 1:return Interpolate_Face_X_Transformed_Weights(block,face_velocities,Transformed(block,X));
      case 2:return Interpolate_Face_Y_Transformed_Weights(block,face_velocities,Transformed(block,X));
      default:assert(axis==3);return Interpolate_Face_Z_Transformed_Weights(block,face_velocities,Transformed(block,X));}}

    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static TV Interpolate_Face_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX)
    {return TV(Interpolate_Face_X_Transformed(block,face_velocities,DX),Interpolate_Face_Y_Transformed(block,face_velocities,DX),
        Interpolate_Face_Z_Transformed(block,face_velocities,DX));}

    template<class T_FACE_LOOKUP>
    static VECTOR<TV,2> Extrema_Face(const T_BLOCK& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const TV& X)
    {return Extrema_Face_Transformed(block,u_min,u_max,Transformed(block,X));}

    template<class T_FACE_LOOKUP>
    static VECTOR<T,2> Extrema_Face_Component(const int axis,const T_BLOCK& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const TV& X)
    {switch(axis){
      case 1:return Extrema_Face_X_Transformed(block,u_min,u_max,Transformed(block,X));
      case 2:return Extrema_Face_Y_Transformed(block,u_min,u_max,Transformed(block,X));
      default:assert(axis==3);return Extrema_Face_Z_Transformed(block,u_min,u_max,Transformed(block,X));}}

    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static VECTOR<TV,2> Extrema_Face_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const TV& DX)
    {VECTOR<T,2> x_extrema=Extrema_Face_X_Transformed(block,u_min,u_max,DX),y_extrema=Extrema_Face_Y_Transformed(block,u_min,u_max,DX),
         z_extrema=Extrema_Face_Z_Transformed(block,u_min,u_max,DX);
    return VECTOR<TV,2>(TV(x_extrema.x,y_extrema.x,z_extrema.x),TV(x_extrema.y,y_extrema.y,z_extrema.y));}

private:
    static TV Transformed(const T_BLOCK& block,const TV& X)
    {return (X-block.Minimum_Corner())*block.One_Over_DX();}

public:
    template<class T_ARRAYS_2>
    static typename T_ARRAYS_2::ELEMENT Interpolate_Cell(const T_BLOCK& block,const T_ARRAYS_2& cell_value,const TV& X)
    {return LINEAR_INTERPOLATION<T,typename T_ARRAYS_2::ELEMENT>::Trilinear(cell_value(block.Cell(0)),cell_value(block.Cell(1)),cell_value(block.Cell(2)),cell_value(block.Cell(3)),
        cell_value(block.Cell(4)),cell_value(block.Cell(5)),cell_value(block.Cell(6)),cell_value(block.Cell(7)),Transformed(block,X));}

//#####################################################################
    static void Block_Transfer(const T_BLOCK& source_block,const T_FACE_ARRAYS_BOOL& source_values,const BLOCK_UNIFORM<GRID<TV> >& destination_block,ARRAY<T,FACE_INDEX<3> >& destination_values);
    TV Interpolate_Face(const TV& X) const;
    // assumes face_velocities are 0 where not valid
    static TV Interpolate_Face_Normalized(const T_BLOCK& block,const T_FACE_ARRAYS& face_velocities,const T_FACE_ARRAYS_BOOL& face_velocities_valid,const TV& X,const TV& default_value=TV());

    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static T Interpolate_Face_X_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX); // between 0 and 1 in the dual cell
    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static T Interpolate_Face_Y_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX); // between 0 and 1 in the dual cell
    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static T Interpolate_Face_Z_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX); // between 0 and 1 in the dual cell
    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static ARRAY<PAIR<FACE_INDEX<3>,T> > Interpolate_Face_X_Transformed_Weights(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX); // between 0 and 1 in the dual cell
    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static ARRAY<PAIR<FACE_INDEX<3>,T> > Interpolate_Face_Y_Transformed_Weights(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX); // between 0 and 1 in the dual cell
    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static ARRAY<PAIR<FACE_INDEX<3>,T> > Interpolate_Face_Z_Transformed_Weights(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX); // between 0 and 1 in the dual cell

    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static VECTOR<T,2> Extrema_Face_X_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const TV& DX); // between 0 and 1 in the dual cell
    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static VECTOR<T,2> Extrema_Face_Y_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const TV& DX); // between 0 and 1 in the dual cell
    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static VECTOR<T,2> Extrema_Face_Z_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const TV& DX); // between 0 and 1 in the dual cell
//#####################################################################
};
}
#endif
