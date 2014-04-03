#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
CUBIC_MN_INTERPOLATION_UNIFORM()
{
    Set_Sharpness();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
~CUBIC_MN_INTERPOLATION_UNIFORM()
{
}
namespace{
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP,class T> T2
From_Base_Node(const CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index)
{
    return cub.cubic_mn_interpolation.Cubic_MN(u(index.x),u(index.x+1),u(index.x+2),u(index.x+3),(X.x-grid.X(index.x+1).x)*grid.one_over_dX.x);
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP,class T> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> >
From_Base_Node_Weights(const CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index)
{
    typedef VECTOR<int,1> TV_INT;
    ARRAY<T> local_weights;
    local_weights=cub.cubic_mn_interpolation.Cubic_MN_Weights((X.x-grid.X(index.x+1).x)*grid.one_over_dX.x);
    ARRAY<PAIR<TV_INT,T> > weights;
    for(int i=0;i<4;i++) weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+i),local_weights(i+1)));
    return weights;
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP,class T> T2
From_Base_Node(const CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index)
{
    VECTOR<T,2> alpha=X-grid.X(index+VECTOR<int,2>::All_Ones_Vector());alpha*=grid.one_over_dX;
    T2 interpolated_in_x[4],value[4];
    for(int b=0;b<4;b++){
        for(int a=0;a<4;a++) value[a]=u(index.x+a,index.y+b);
        interpolated_in_x[b]=cub.cubic_mn_interpolation.Cubic_MN(value[0],value[1],value[2],value[3],alpha.x);}
    return cub.cubic_mn_interpolation.Cubic_MN(interpolated_in_x[0],interpolated_in_x[1],interpolated_in_x[2],interpolated_in_x[3],alpha.y);
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP,class T> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> >
From_Base_Node_Weights(const CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index)
{
    typedef VECTOR<int,2> TV_INT;
    VECTOR<T,2> alpha=X-grid.X(index+VECTOR<int,2>::All_Ones_Vector());alpha*=grid.one_over_dX;
    ARRAY<PAIR<TV_INT,T> > weights;
    ARRAY<T> local_weights_x=cub.cubic_mn_interpolation.Cubic_MN_Weights(alpha.x);
    ARRAY<T> local_weights_y=cub.cubic_mn_interpolation.Cubic_MN_Weights(alpha.y);
    for(int b=0;b<4;b++) for(int a=0;a<4;a++) weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+a,index.y+b),local_weights_y(b+1)*local_weights_x(a+1)));
    return weights;
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP,class T> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> >
From_Base_Node_Weights(const CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index)
{
    typedef VECTOR<int,3> TV_INT;
    ARRAY<PAIR<TV_INT,T> > weights;
    PHYSBAM_NOT_IMPLEMENTED();
    return weights;
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP,class T> T2
From_Base_Node(const CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,3> >& u,const VECTOR<T,3>& X,
    const VECTOR<int,3>& index)
{
    VECTOR<T,3> alpha=X-grid.X(index+VECTOR<int,3>::All_Ones_Vector());alpha*=grid.one_over_dX;
    T2 interpolated_in_x[4],interpolated_in_y[4],value[4];
    for(int c=0;c<4;c++){
        for(int b=0;b<4;b++){
            for(int a=0;a<4;a++) value[a]=u(index.x+a,index.y+b,index.z+c);
            interpolated_in_x[b]=cub.cubic_mn_interpolation.Cubic_MN(value[0],value[1],value[2],value[3],alpha.x);}
        interpolated_in_y[c]=cub.cubic_mn_interpolation.Cubic_MN(interpolated_in_x[0],interpolated_in_x[1],interpolated_in_x[2],interpolated_in_x[3],alpha.y);}
    return cub.cubic_mn_interpolation.Cubic_MN(interpolated_in_y[0],interpolated_in_y[1],interpolated_in_y[2],interpolated_in_y[3],alpha.z);
}
}
//#####################################################################
// Function Clamped_To_Array
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE
{
    return From_Base_Node(*this,grid,u,X,INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::Clamped_Index_Interior_End_Minus_One(grid,u,X)-TV_INT::All_Ones_Vector());
}
//#####################################################################
// Function Clamped_To_Array_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> > CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Clamped_To_Array_Weights(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE
{
    return From_Base_Node_Weights(*this,grid,u,X,INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::Clamped_Index_Interior_End_Minus_One(grid,u,X)-TV_INT::All_Ones_Vector());
}
namespace{
//#####################################################################
// Function From_Base_Node_Periodic
//#####################################################################
template<class T,class T_GRID,class T2,class T_FACE_LOOKUP> T2
    From_Base_Node_Periodic(const CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u,const VECTOR<T,2>& X,
        const VECTOR<int,2>& index)
{
    VECTOR<T,2> alpha=X-grid.X(index+VECTOR<int,2>::All_Ones_Vector());alpha*=grid.one_over_dX;
    T2 interpolated_in_x[4],value[4];
    for(int b=0;b<4;b++){
        for(int a=0;a<4;a++){
            int x=index.x+a-1,y=index.y+b-1;
            while(x<0) x+=grid.counts.x;
            if(x>=grid.counts.x)x=(x&(grid.counts.x-1));
            while(y<0) y+=grid.counts.y;
            if(y>=grid.counts.y)y=(y&(grid.counts.y-1));
            value[a]=u(x+1,y+1);}
        interpolated_in_x[b]=cub.cubic_mn_interpolation.Cubic_MN(value[0],value[1],value[2],value[3],alpha.x);}
    return cub.cubic_mn_interpolation.Cubic_MN(interpolated_in_x[0],interpolated_in_x[1],interpolated_in_x[2],interpolated_in_x[3],alpha.y);
}
template<class T,class TV,class T_GRID,class T2,class T_FACE_LOOKUP,int d> T2
    From_Base_Node_Periodic(const CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>& cub,const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,d> >& u,const VECTOR<T,d>& X,
        const VECTOR<int,d>& index)
{
    PHYSBAM_FATAL_ERROR();
}
}
//#####################################################################
// Function Periodic
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 CUBIC_MN_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Periodic(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const
{
    return From_Base_Node_Periodic(*this,grid,u,X,grid.Cell(X,0)-TV_INT::All_Ones_Vector());
}
template class CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >;
template class CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >;
template class CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >;
template class CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >;
template class CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >;
template class CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >;
template class CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >;
#endif
