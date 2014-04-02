//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_2D_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_2D_HELPER_DEFINITIONS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
using namespace PhysBAM;
template<class T_GRID> LINEAR_INTERPOLATION_MAC_2D_HELPER<T_GRID>::
LINEAR_INTERPOLATION_MAC_2D_HELPER(const T_BLOCK& block,const T_FACE_ARRAYS_SCALAR& face_velocities)
    :base(block.Minimum_Corner()),center(block.Center()),one_over_DX(block.One_Over_DX())
{
    typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP face_velocities_lookup(face_velocities);
    static const int rotated_face_x[6]={0,1,2,3,4,5};
    u2=block.Face_X_Value(face_velocities_lookup,rotated_face_x[1]);u5=block.Face_X_Value(face_velocities_lookup,rotated_face_x[4]);
    slope_u12=one_over_DX.x*(u2-block.Face_X_Value(face_velocities_lookup,rotated_face_x[0]));
    slope_u23=one_over_DX.x*(block.Face_X_Value(face_velocities_lookup,rotated_face_x[2])-u2);
    slope_u45=one_over_DX.x*(u5-block.Face_X_Value(face_velocities_lookup,rotated_face_x[3]));
    slope_u56=one_over_DX.x*(block.Face_X_Value(face_velocities_lookup,rotated_face_x[5])-u5);
    static const int rotated_face_y[6]={0,2,4,1,3,5};
    v2=block.Face_Y_Value(face_velocities_lookup,rotated_face_y[1]);v5=block.Face_Y_Value(face_velocities_lookup,rotated_face_y[4]);
    slope_v12=one_over_DX.y*(v2-block.Face_Y_Value(face_velocities_lookup,rotated_face_y[0]));
    slope_v23=one_over_DX.y*(block.Face_Y_Value(face_velocities_lookup,rotated_face_y[2])-v2);
    slope_v45=one_over_DX.y*(v5-block.Face_Y_Value(face_velocities_lookup,rotated_face_y[3]));
    slope_v56=one_over_DX.y*(block.Face_Y_Value(face_velocities_lookup,rotated_face_y[5])-v5);
}
template<class T_GRID> VECTOR<typename T_GRID::VECTOR_T::SCALAR,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<T_GRID>::
Interpolate_Face(const VECTOR<T,2>& X) const
{
    VECTOR<T,2> yx(X.y,X.x);
    return VECTOR<T,2>(X.x<center.x?LINEAR_INTERPOLATION<T,T>::Bilinear(u2,u5,one_over_DX.y,center.x,base.y,slope_u12,slope_u45,X)
                                    :LINEAR_INTERPOLATION<T,T>::Bilinear(u2,u5,one_over_DX.y,center.x,base.y,slope_u23,slope_u56,X),
                        X.y<center.y?LINEAR_INTERPOLATION<T,T>::Bilinear(v2,v5,one_over_DX.x,center.y,base.x,slope_v12,slope_v45,yx)
                                    :LINEAR_INTERPOLATION<T,T>::Bilinear(v2,v5,one_over_DX.x,center.y,base.x,slope_v23,slope_v56,yx));
}
// assumes face_velocities are 0 where not valid
template<class T_GRID> VECTOR<typename T_GRID::VECTOR_T::SCALAR,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<T_GRID>::
Interpolate_Face_Normalized(const T_BLOCK& block,const T_FACE_ARRAYS_SCALAR& face_velocities,const T_FACE_ARRAYS_BOOL& face_velocities_valid,const VECTOR<T,2>& X,const VECTOR<T,2>& default_value)
{
    static const GRID<TV> valid_values_grid=GRID<TV>(2,2,0,1,0,1).Get_MAC_Grid_At_Regular_Positions();
    static const BLOCK_UNIFORM<GRID<TV> > valid_values_block(valid_values_grid,VECTOR<int,2>(2,2));
    ARRAY<T,FACE_INDEX<2> > valid_values(valid_values_grid);Block_Transfer(block,face_velocities_valid,valid_values_block,valid_values);
    VECTOR<T,2> DX=Transformed(block,X),velocity=Interpolate_Face_Transformed(block,face_velocities,DX),weight=Interpolate_Face_Transformed(valid_values_block,valid_values,DX);
    return VECTOR<T,2>(weight.x==0?default_value.x:velocity.x/weight.x,weight.y==0?default_value.y:velocity.y/weight.y);
}
template class LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >;
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > const&,VECTOR<float,2> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > const&,VECTOR<float,2> const&);
template float LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > const&,VECTOR<float,2> const&);
template float LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > const&,VECTOR<float,2> const&);
template ARRAY<PAIR<FACE_INDEX<2>,float> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > const&,VECTOR<float,2> const&);
template ARRAY<PAIR<FACE_INDEX<2>,float> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > const&,VECTOR<float,2> const&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > const&,VECTOR<double,2> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > const&,VECTOR<double,2> const&);
template class LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >;
template double LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > const&,VECTOR<double,2> const&);
template double LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > const&,VECTOR<double,2> const&);
template ARRAY<PAIR<FACE_INDEX<2>,double> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > const&,VECTOR<double,2> const&);
template ARRAY<PAIR<FACE_INDEX<2>,double> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > const&,VECTOR<double,2> const&);
#endif
