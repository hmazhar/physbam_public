//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER_DEFINITIONS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
using namespace PhysBAM;
template<class T_GRID> LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
LINEAR_INTERPOLATION_MAC_3D_HELPER(const T_BLOCK& block,const T_FACE_ARRAYS& face_velocities)
    :base(block.Minimum_Corner()),center(block.Center()),one_over_DX(block.One_Over_DX())
{
    typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP face_velocities_lookup(face_velocities);
    static const int rotated_face_x[12]={0,1,2,3,4,5,6,7,8,9,10,11};
    u2=block.Face_X_Value(face_velocities_lookup,rotated_face_x[1]);u5=block.Face_X_Value(face_velocities_lookup,rotated_face_x[4]);
    u8=block.Face_X_Value(face_velocities_lookup,rotated_face_x[7]);u11=block.Face_X_Value(face_velocities_lookup,rotated_face_x[10]);
    slope_u12=one_over_DX.x*(u2-block.Face_X_Value(face_velocities_lookup,rotated_face_x[0]));
    slope_u23=one_over_DX.x*(block.Face_X_Value(face_velocities_lookup,rotated_face_x[2])-u2);
    slope_u45=one_over_DX.x*(u5-block.Face_X_Value(face_velocities_lookup,rotated_face_x[3]));
    slope_u56=one_over_DX.x*(block.Face_X_Value(face_velocities_lookup,rotated_face_x[5])-u5);
    slope_u78=one_over_DX.x*(u8-block.Face_X_Value(face_velocities_lookup,rotated_face_x[6]));
    slope_u89=one_over_DX.x*(block.Face_X_Value(face_velocities_lookup,rotated_face_x[8])-u8);
    slope_u10_11=one_over_DX.x*(u11-block.Face_X_Value(face_velocities_lookup,rotated_face_x[9]));
    slope_u11_12=one_over_DX.x*(block.Face_X_Value(face_velocities_lookup,rotated_face_x[11])-u11);
    static const int rotated_face_y[12]={0,2,4,1,3,5,6,8,10,7,9,11};
    v2=block.Face_Y_Value(face_velocities_lookup,rotated_face_y[1]);v5=block.Face_Y_Value(face_velocities_lookup,rotated_face_y[4]);
    v8=block.Face_Y_Value(face_velocities_lookup,rotated_face_y[7]);v11=block.Face_Y_Value(face_velocities_lookup,rotated_face_y[10]);
    slope_v12=one_over_DX.y*(v2-block.Face_Y_Value(face_velocities_lookup,rotated_face_y[0]));
    slope_v23=one_over_DX.y*(block.Face_Y_Value(face_velocities_lookup,rotated_face_y[2])-v2);
    slope_v45=one_over_DX.y*(v5-block.Face_Y_Value(face_velocities_lookup,rotated_face_y[3]));
    slope_v56=one_over_DX.y*(block.Face_Y_Value(face_velocities_lookup,rotated_face_y[5])-v5);
    slope_v78=one_over_DX.y*(v8-block.Face_Y_Value(face_velocities_lookup,rotated_face_y[6]));
    slope_v89=one_over_DX.y*(block.Face_Y_Value(face_velocities_lookup,rotated_face_y[8])-v8);
    slope_v10_11=one_over_DX.y*(v11-block.Face_Y_Value(face_velocities_lookup,rotated_face_y[9]));
    slope_v11_12=one_over_DX.y*(block.Face_Y_Value(face_velocities_lookup,rotated_face_y[11])-v11);
    static const int rotated_face_z[12]={0,4,8,1,5,9,2,6,10,3,7,11};
    w2=block.Face_Z_Value(face_velocities_lookup,rotated_face_z[1]);w5=block.Face_Z_Value(face_velocities_lookup,rotated_face_z[4]);
    w8=block.Face_Z_Value(face_velocities_lookup,rotated_face_z[7]);w11=block.Face_Z_Value(face_velocities_lookup,rotated_face_z[10]);
    slope_w12=one_over_DX.z*(w2-block.Face_Z_Value(face_velocities_lookup,rotated_face_z[0]));
    slope_w23=one_over_DX.z*(block.Face_Z_Value(face_velocities_lookup,rotated_face_z[2])-w2);
    slope_w45=one_over_DX.z*(w5-block.Face_Z_Value(face_velocities_lookup,rotated_face_z[3]));
    slope_w56=one_over_DX.z*(block.Face_Z_Value(face_velocities_lookup,rotated_face_z[5])-w5);
    slope_w78=one_over_DX.z*(w8-block.Face_Z_Value(face_velocities_lookup,rotated_face_z[6]));
    slope_w89=one_over_DX.z*(block.Face_Z_Value(face_velocities_lookup,rotated_face_z[8])-w8);
    slope_w10_11=one_over_DX.z*(w11-block.Face_Z_Value(face_velocities_lookup,rotated_face_z[9]));
    slope_w11_12=one_over_DX.z*(block.Face_Z_Value(face_velocities_lookup,rotated_face_z[11])-w11);
}
template<class T_GRID> LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
~LINEAR_INTERPOLATION_MAC_3D_HELPER()
{
}
template<class T_GRID> typename T_GRID::VECTOR_T LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
Interpolate_Face(const TV& X) const
{
    TV yxz(X.y,X.x,X.z),zxy(X.z,X.x,X.y);
    return TV(X.x<center.x?LINEAR_INTERPOLATION<T,T>::Trilinear(u2,u5,u8,u11,one_over_DX.y,one_over_DX.z,center.x,base.y,base.z,slope_u12,slope_u45,slope_u78,slope_u10_11,X)
        :LINEAR_INTERPOLATION<T,T>::Trilinear(u2,u5,u8,u11,one_over_DX.y,one_over_DX.z,center.x,base.y,base.z,slope_u23,slope_u56,slope_u89,slope_u11_12,X),
        X.y<center.y?LINEAR_INTERPOLATION<T,T>::Trilinear(v2,v5,v8,v11,one_over_DX.x,one_over_DX.z,center.y,base.x,base.z,slope_v12,slope_v45,slope_v78,slope_v10_11,yxz)
        :LINEAR_INTERPOLATION<T,T>::Trilinear(v2,v5,v8,v11,one_over_DX.x,one_over_DX.z,center.y,base.x,base.z,slope_v23,slope_v56,slope_v89,slope_v11_12,yxz),
        X.z<center.z?LINEAR_INTERPOLATION<T,T>::Trilinear(w2,w5,w8,w11,one_over_DX.x,one_over_DX.y,center.z,base.x,base.y,slope_w12,slope_w45,slope_w78,slope_w10_11,zxy)
        :LINEAR_INTERPOLATION<T,T>::Trilinear(w2,w5,w8,w11,one_over_DX.x,one_over_DX.y,center.z,base.x,base.y,slope_w23,slope_w56,slope_w89,slope_w11_12,zxy));
}
template<class T_GRID> void LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
Block_Transfer(const T_BLOCK& source_block,const T_FACE_ARRAYS_BOOL& source_values,const BLOCK_UNIFORM<GRID<TV> >& destination_block,ARRAY<T,FACE_INDEX<3> >& destination_values)
{
    for(int i=0;i<T_GRID::number_of_faces_per_block/T_GRID::dimension;i++){
        destination_block.Face_X_Reference(destination_values,i)=(T)source_block.Face_X_Value(source_values,i);
        destination_block.Face_Y_Reference(destination_values,i)=(T)source_block.Face_Y_Value(source_values,i);
        destination_block.Face_Z_Reference(destination_values,i)=(T)source_block.Face_Z_Value(source_values,i);}
}
template<class T_GRID> typename T_GRID::VECTOR_T LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
Interpolate_Face_Normalized(const T_BLOCK& block,const T_FACE_ARRAYS& face_velocities,const T_FACE_ARRAYS_BOOL& face_velocities_valid,const TV& X,const TV& default_value)
{
    static const GRID<TV> valid_values_grid=GRID<TV>(2,2,2,0,1,0,1,0,1).Get_MAC_Grid_At_Regular_Positions();
    static const BLOCK_UNIFORM<GRID<TV> > valid_values_block(valid_values_grid,VECTOR<int,3>(2,2,2));
    ARRAY<T,FACE_INDEX<3> > valid_values(valid_values_grid);Block_Transfer(block,face_velocities_valid,valid_values_block,valid_values);
    TV DX=Transformed(block,X),velocity=Interpolate_Face_Transformed(block,face_velocities,DX),weight=Interpolate_Face_Transformed(valid_values_block,valid_values,DX);
    return TV(weight.x==0?default_value.x:velocity.x/weight.x,weight.y==0?default_value.y:velocity.y/weight.y,weight.z==0?default_value.z:velocity.z/weight.z);
}
template class LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >;
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > const&,VECTOR<float,3> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > const&,VECTOR<float,3> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Extrema_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Z_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > const&,VECTOR<float,3> const&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >;
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > const&,VECTOR<double,3> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > const&,VECTOR<double,3> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Extrema_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Z_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > const&,VECTOR<double,3> const&);
#endif
