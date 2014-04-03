//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER.h>
#include <PhysBAM_Tools/Math_Tools/Componentwise_Min_Max.h>
using namespace PhysBAM;
template<class T_GRID> template<class T_BLOCK_2,class T_FACE_LOOKUP> typename T_GRID::VECTOR_T::SCALAR LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
Interpolate_Face_X_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX) // between 0 and 1 in the dual cell
{
    if(DX.x<.5) return LINEAR_INTERPOLATION<T,T>::Trilinear(
        block.Face_X_Value(face_velocities,0),block.Face_X_Value(face_velocities,1),
        block.Face_X_Value(face_velocities,3),block.Face_X_Value(face_velocities,4),
        block.Face_X_Value(face_velocities,6),block.Face_X_Value(face_velocities,7),
        block.Face_X_Value(face_velocities,9),block.Face_X_Value(face_velocities,10),TV(DX.x+(T).5,DX.y,DX.z));
    else return LINEAR_INTERPOLATION<T,T>::Trilinear(
        block.Face_X_Value(face_velocities,1),block.Face_X_Value(face_velocities,2),
        block.Face_X_Value(face_velocities,4),block.Face_X_Value(face_velocities,5),
        block.Face_X_Value(face_velocities,7),block.Face_X_Value(face_velocities,8),
        block.Face_X_Value(face_velocities,10),block.Face_X_Value(face_velocities,11),TV(DX.x-(T).5,DX.y,DX.z));
}
template<class T_GRID> template<class T_BLOCK_2,class T_FACE_LOOKUP> typename T_GRID::VECTOR_T::SCALAR LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
Interpolate_Face_Y_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX) // between 0 and 1 in the dual cell
{
    if(DX.y<.5) return LINEAR_INTERPOLATION<T,T>::Trilinear(
        block.Face_Y_Value(face_velocities,0),block.Face_Y_Value(face_velocities,1),
        block.Face_Y_Value(face_velocities,2),block.Face_Y_Value(face_velocities,3),
        block.Face_Y_Value(face_velocities,6),block.Face_Y_Value(face_velocities,7),
        block.Face_Y_Value(face_velocities,8),block.Face_Y_Value(face_velocities,9),TV(DX.x,DX.y+(T).5,DX.z));
    else return LINEAR_INTERPOLATION<T,T>::Trilinear(
        block.Face_Y_Value(face_velocities,2),block.Face_Y_Value(face_velocities,3),
        block.Face_Y_Value(face_velocities,4),block.Face_Y_Value(face_velocities,5),
        block.Face_Y_Value(face_velocities,8),block.Face_Y_Value(face_velocities,9),
        block.Face_Y_Value(face_velocities,10),block.Face_Y_Value(face_velocities,11),TV(DX.x,DX.y-(T).5,DX.z));
}
template<class T_GRID> template<class T_BLOCK_2,class T_FACE_LOOKUP> typename T_GRID::VECTOR_T::SCALAR LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
Interpolate_Face_Z_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX) // between 0 and 1 in the dual cell
{
    if(DX.z<.5) return LINEAR_INTERPOLATION<T,T>::Trilinear(
        block.Face_Z_Value(face_velocities,0),block.Face_Z_Value(face_velocities,1),
        block.Face_Z_Value(face_velocities,2),block.Face_Z_Value(face_velocities,3),
        block.Face_Z_Value(face_velocities,4),block.Face_Z_Value(face_velocities,5),
        block.Face_Z_Value(face_velocities,6),block.Face_Z_Value(face_velocities,7),TV(DX.x,DX.y,DX.z+(T).5));
    else return LINEAR_INTERPOLATION<T,T>::Trilinear(
        block.Face_Z_Value(face_velocities,4),block.Face_Z_Value(face_velocities,5),
        block.Face_Z_Value(face_velocities,6),block.Face_Z_Value(face_velocities,7),
        block.Face_Z_Value(face_velocities,8),block.Face_Z_Value(face_velocities,9),
        block.Face_Z_Value(face_velocities,10),block.Face_Z_Value(face_velocities,11),TV(DX.x,DX.y,DX.z-(T).5));
}
template<class T_GRID> template<class T_BLOCK_2,class T_FACE_LOOKUP> ARRAY<PAIR<FACE_INDEX<3>,typename T_GRID::VECTOR_T::SCALAR> > LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
Interpolate_Face_X_Transformed_Weights(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX) // between 0 and 1 in the dual cell
{
    ARRAY<PAIR<FACE_INDEX<3>,T> > weights;
    if(DX.x<.5){VECTOR<T,3> w(DX.x+(T).5,DX.y,DX.z);
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(0),(1-w.x)*(1-w.y)*(1-w.z)));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(1),w.x*(1-w.y)*(1-w.z)));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(3),(1-w.x)*w.y*(1-w.z)));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(4),w.x*w.y*(1-w.z)));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(6),(1-w.x)*(1-w.y)*w.z));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(7),w.x*(1-w.y)*w.z));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(9),(1-w.x)*w.y*w.z));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(10),w.x*w.y*w.z));return weights;}
    else{VECTOR<T,3> w(DX.x-(T).5,DX.y,DX.z);
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(1),(1-w.x)*(1-w.y)*(1-w.z)));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(2),w.x*(1-w.y)*(1-w.z)));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(4),(1-w.x)*w.y*(1-w.z)));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(5),w.x*w.y*(1-w.z)));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(7),(1-w.x)*(1-w.y)*w.z));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(8),w.x*(1-w.y)*w.z));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(10),(1-w.x)*w.y*w.z));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_X_Index(11),w.x*w.y*w.z));return weights;}
}
template<class T_GRID> template<class T_BLOCK_2,class T_FACE_LOOKUP> ARRAY<PAIR<FACE_INDEX<3>,typename T_GRID::VECTOR_T::SCALAR> > LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
Interpolate_Face_Y_Transformed_Weights(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX) // between 0 and 1 in the dual cell
{
    ARRAY<PAIR<FACE_INDEX<3>,T> > weights;
    if(DX.y<.5){VECTOR<T,3> w(DX.x,DX.y+(T).5,DX.z);
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(0),(1-w.x)*(1-w.y)*(1-w.z)));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(1),w.x*(1-w.y)*(1-w.z)));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(2),(1-w.x)*w.y*(1-w.z)));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(3),w.x*w.y*(1-w.z)));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(6),(1-w.x)*(1-w.y)*w.z));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(7),w.x*(1-w.y)*w.z));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(8),(1-w.x)*w.y*w.z));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(9),w.x*w.y*w.z));return weights;}
    else{VECTOR<T,3> w(DX.x,DX.y-(T).5,DX.z);
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(2),(1-w.x)*(1-w.y)*(1-w.z)));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(3),w.x*(1-w.y)*(1-w.z)));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(4),(1-w.x)*w.y*(1-w.z)));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(5),w.x*w.y*(1-w.z)));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(8),(1-w.x)*(1-w.y)*w.z));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(9),w.x*(1-w.y)*w.z));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(10),(1-w.x)*w.y*w.z));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Y_Index(11),w.x*w.y*w.z));return weights;}
}
template<class T_GRID> template<class T_BLOCK_2,class T_FACE_LOOKUP> ARRAY<PAIR<FACE_INDEX<3>,typename T_GRID::VECTOR_T::SCALAR> > LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
Interpolate_Face_Z_Transformed_Weights(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const TV& DX) // between 0 and 1 in the dual cell
{ 
    ARRAY<PAIR<FACE_INDEX<3>,T> > weights;
    if(DX.z<.5){VECTOR<T,3> w(DX.x,DX.y,DX.z+(T).5);
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(0),(1-w.x)*(1-w.y)*(1-w.z)));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(1),w.x*(1-w.y)*(1-w.z)));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(2),(1-w.x)*w.y*(1-w.z)));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(3),w.x*w.y*(1-w.z)));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(4),(1-w.x)*(1-w.y)*w.z));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(5),w.x*(1-w.y)*w.z));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(6),(1-w.x)*w.y*w.z));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(7),w.x*w.y*w.z));return weights;}
    else{VECTOR<T,3> w(DX.x,DX.y,DX.z-(T).5);
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(4),(1-w.x)*(1-w.y)*(1-w.z)));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(5),w.x*(1-w.y)*(1-w.z)));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(6),(1-w.x)*w.y*(1-w.z)));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(7),w.x*w.y*(1-w.z)));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(8),(1-w.x)*(1-w.y)*w.z));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(9),w.x*(1-w.y)*w.z));
        weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(10),(1-w.x)*w.y*w.z));weights.Append(PAIR<FACE_INDEX<3>,T>(block.Face_Z_Index(11),w.x*w.y*w.z));return weights;}
}
template<class T_GRID> template<class T_BLOCK_2,class T_FACE_LOOKUP> VECTOR<typename T_GRID::VECTOR_T::SCALAR,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
Extrema_Face_X_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const TV& DX) // between 0 and 1 in the dual cell
{
    if(DX.x<.5) return VECTOR<T,2>(
         Componentwise_Min(block.Face_X_Value(u_min,0),block.Face_X_Value(u_min,1),block.Face_X_Value(u_min,3),block.Face_X_Value(u_min,4),
             block.Face_X_Value(u_min,6),block.Face_X_Value(u_min,7),block.Face_X_Value(u_min,9),block.Face_X_Value(u_min,10)),
         Componentwise_Max(block.Face_X_Value(u_max,0),block.Face_X_Value(u_max,1),block.Face_X_Value(u_max,3),block.Face_X_Value(u_max,4),
             block.Face_X_Value(u_max,6),block.Face_X_Value(u_max,7),block.Face_X_Value(u_max,9),block.Face_X_Value(u_max,10)));
    else return VECTOR<T,2>(
         Componentwise_Min(block.Face_X_Value(u_min,1),block.Face_X_Value(u_min,2),block.Face_X_Value(u_min,4),block.Face_X_Value(u_min,5),
             block.Face_X_Value(u_min,7),block.Face_X_Value(u_min,8),block.Face_X_Value(u_min,10),block.Face_X_Value(u_min,11)),
         Componentwise_Max(block.Face_X_Value(u_max,1),block.Face_X_Value(u_max,2),block.Face_X_Value(u_max,4),block.Face_X_Value(u_max,5),
             block.Face_X_Value(u_max,7),block.Face_X_Value(u_max,8),block.Face_X_Value(u_max,10),block.Face_X_Value(u_max,11)));
}
template<class T_GRID> template<class T_BLOCK_2,class T_FACE_LOOKUP> VECTOR<typename T_GRID::VECTOR_T::SCALAR,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
Extrema_Face_Y_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const TV& DX) // between 0 and 1 in the dual cell
{
    if(DX.y<.5) return VECTOR<T,2>(
        Componentwise_Min(block.Face_Y_Value(u_min,0),block.Face_Y_Value(u_min,1),block.Face_Y_Value(u_min,2),block.Face_Y_Value(u_min,3),
            block.Face_Y_Value(u_min,6),block.Face_Y_Value(u_min,7),block.Face_Y_Value(u_min,8),block.Face_Y_Value(u_min,9)),
        Componentwise_Max(block.Face_Y_Value(u_max,0),block.Face_Y_Value(u_max,1),block.Face_Y_Value(u_max,2),block.Face_Y_Value(u_max,3),
            block.Face_Y_Value(u_max,6),block.Face_Y_Value(u_max,7),block.Face_Y_Value(u_max,8),block.Face_Y_Value(u_max,9)));
    else return VECTOR<T,2>(
        Componentwise_Min(block.Face_Y_Value(u_min,2),block.Face_Y_Value(u_min,3),block.Face_Y_Value(u_min,4),block.Face_Y_Value(u_min,5),
            block.Face_Y_Value(u_min,8),block.Face_Y_Value(u_min,9),block.Face_Y_Value(u_min,10),block.Face_Y_Value(u_min,11)),
        Componentwise_Max(block.Face_Y_Value(u_max,2),block.Face_Y_Value(u_max,3),block.Face_Y_Value(u_max,4),block.Face_Y_Value(u_max,5),
            block.Face_Y_Value(u_max,8),block.Face_Y_Value(u_max,9),block.Face_Y_Value(u_max,10),block.Face_Y_Value(u_max,11)));
}
template<class T_GRID> template<class T_BLOCK_2,class T_FACE_LOOKUP> VECTOR<typename T_GRID::VECTOR_T::SCALAR,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<T_GRID>::
Extrema_Face_Z_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const TV& DX) // between 0 and 1 in the dual cell
{
    if(DX.z<.5) return VECTOR<T,2>(
        Componentwise_Min(block.Face_Z_Value(u_min,0),block.Face_Z_Value(u_min,1),block.Face_Z_Value(u_min,2),block.Face_Z_Value(u_min,3),
            block.Face_Z_Value(u_min,4),block.Face_Z_Value(u_min,5),block.Face_Z_Value(u_min,6),block.Face_Z_Value(u_min,7)),
        Componentwise_Max(block.Face_Z_Value(u_max,0),block.Face_Z_Value(u_max,1),block.Face_Z_Value(u_max,2),block.Face_Z_Value(u_max,3),
            block.Face_Z_Value(u_max,4),block.Face_Z_Value(u_max,5),block.Face_Z_Value(u_max,6),block.Face_Z_Value(u_max,7)));
    else return VECTOR<T,2>(
        Componentwise_Min(block.Face_Z_Value(u_min,4),block.Face_Z_Value(u_min,5),block.Face_Z_Value(u_min,6),block.Face_Z_Value(u_min,7),
            block.Face_Z_Value(u_min,8),block.Face_Z_Value(u_min,9),block.Face_Z_Value(u_min,10),block.Face_Z_Value(u_min,11)),
        Componentwise_Max(block.Face_Z_Value(u_max,4),block.Face_Z_Value(u_max,5),block.Face_Z_Value(u_max,6),block.Face_Z_Value(u_max,7),
            block.Face_Z_Value(u_max,8),block.Face_Z_Value(u_max,9),block.Face_Z_Value(u_max,10),block.Face_Z_Value(u_max,11)));
}
