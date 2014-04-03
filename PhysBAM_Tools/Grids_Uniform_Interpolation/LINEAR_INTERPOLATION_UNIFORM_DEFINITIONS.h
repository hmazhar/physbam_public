#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Function Clamped_To_Array_No_Extrema
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Clamped_To_Array_No_Extrema(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const
{
    TV_INT index=INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::Clamped_Index_End_Minus_One(grid,u,X);
    TV X_clamped=clamp(X,grid.X(index),grid.X(index+1));
    return From_Base_Node(grid,u,X_clamped,index);
}
//#####################################################################
// Function Clamped_To_Array
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const
{
    TV_INT index=INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::Clamped_Index_End_Minus_One(grid,u,X);
    return From_Base_Node(grid,u,X,index);
}
//#####################################################################
// Function Clamped_To_Array_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> > LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Clamped_To_Array_Weights(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const
{
    TV_INT index=INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::Clamped_Index_End_Minus_One(grid,u,X);
    return From_Base_Node_Weights(grid,u,X,index);
}
namespace{
//#####################################################################
// Function Compute_Inverse_Map_Helper
//#####################################################################
template<class T>
void Compute_Inverse_Map_Helper(const GRID<VECTOR<T,1> >& domain_grid,const ARRAY<T,VECTOR<int,1> >& function,const GRID<VECTOR<T,1> >& range_grid,ARRAY<T,VECTOR<int,1> >& inverse_function)
{
    int domain_i=1;
    T xmin=domain_grid.Axis_X(1,1),xmax=domain_grid.Axis_X(domain_grid.counts.x,1);
    for(int i=1;i<=range_grid.counts.x;i++){
        T function_value=range_grid.Axis_X(i,1);
        while(domain_i<domain_grid.counts.x-1 && function(domain_i+1)<function_value) domain_i++;
        inverse_function(i)=clamp(domain_grid.Axis_X(domain_i,1)+domain_grid.dX.x*(function_value-function(domain_i))/(function(domain_i+1)-function(domain_i)),xmin,xmax);}
}
//#####################################################################
// Function Compute_Inverse_Map_Helper
//#####################################################################
template<class T,class T2,int d>
void Compute_Inverse_Map_Helper(const GRID<VECTOR<T,d> >& domain_grid,const ARRAY<T2,VECTOR<int,1> >& function,const GRID<VECTOR<T,d> >& range_grid,ARRAY<T2,VECTOR<int,1> >& inverse_function)
{
    PHYSBAM_FATAL_ERROR();
}
}
//#####################################################################
// Function Compute_Inverse_Map
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> void LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Compute_Inverse_Map(const GRID<TV>& domain_grid,const ARRAY<T2,VECTOR<int,1> >& function,const GRID<TV>& range_grid,ARRAY<T2,VECTOR<int,1> >& inverse_function)
{
    Compute_Inverse_Map_Helper(domain_grid,function,range_grid,inverse_function);
}
//#####################################################################
// Function Extrema_Clamped_To_Array
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> VECTOR<T2,2> LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Extrema_Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u_min,const T_ARRAYS_T2& u_max,const TV& X) const
{
    TV_INT index=INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::Clamped_Index_End_Minus_One(grid,u_min,X);
    return Extrema_From_Base_Node(grid,u_min,u_max,X,index);
}
namespace{
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T,class T2> T2
From_Base_Node_Helper(const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index)
{
    int i=index.x;
    return LINEAR_INTERPOLATION<T,T2>::Linear(u(index),u(i+1),(X-grid.X(index))*grid.one_over_dX);
}
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T,class T2> T2
From_Base_Node_Helper(const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index)
{
    int i=index.x,j=index.y;
    return LINEAR_INTERPOLATION<T,T2>::Bilinear(u(index),u(i+1,j),u(i,j+1),u(i+1,j+1),(X-grid.X(index))*grid.one_over_dX);
}
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T,class T2> T2
From_Base_Node_Helper(const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index)
{
    const T2 *base=&u(index),*basep1=base+u.counts.y*u.counts.z;
    return LINEAR_INTERPOLATION<T,T2>::Trilinear(*base,*basep1,base[u.counts.z],basep1[u.counts.z],base[1],basep1[1],base[u.counts.z+1],basep1[u.counts.z+1],(X-grid.X(index))*grid.one_over_dX);
}
}
namespace{
//#####################################################################
// Function From_Base_Node_Weights_Helper
//#####################################################################
template<class T,class T2> ARRAY<PAIR<VECTOR<int,1>,T> >
From_Base_Node_Weights_Helper(const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index)
{
    typedef VECTOR<int,1> TV_INT;
    int i=index.x;
    T x=((X-grid.X(index))*grid.one_over_dX).x;
    ARRAY<PAIR<TV_INT,T> > weights;
    weights.Append(PAIR<TV_INT,T>(index,1-x));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1),x));
    return weights;
}
//#####################################################################
// Function From_Base_Node_Weights_Helper
//#####################################################################
template<class T,class T2> ARRAY<PAIR<VECTOR<int,2>,T> >
From_Base_Node_Weights_Helper(const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index)
{
    typedef VECTOR<int,2> TV_INT;typedef VECTOR<T,2> TV;
    int i=index.x;
    int j=index.y;
    TV w=(X-grid.X(index))*grid.one_over_dX;
    ARRAY<PAIR<TV_INT,T> > weights;
    weights.Append(PAIR<TV_INT,T>(index,(1-w.x)*(1-w.y)));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1,j),w.x*(1-w.y)));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i,j+1),(1-w.x)*w.y));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1,j+1),w.x*w.y));
    return weights;
}
//#####################################################################
// Function From_Base_Node_Weights_Helper
//#####################################################################
template<class T,class T2> ARRAY<PAIR<VECTOR<int,3>,T> >
From_Base_Node_Weights_Helper(const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index)
{
    typedef VECTOR<int,3> TV_INT;typedef VECTOR<T,3> TV;
    int i=index.x;
    int j=index.y;
    int k=index.z;
    TV w=(X-grid.X(index))*grid.one_over_dX;
    ARRAY<PAIR<TV_INT,T> > weights;
    weights.Append(PAIR<TV_INT,T>(index,(1-w.x)*(1-w.y)*(1-w.z)));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1,j,k),w.x*(1-w.y)*(1-w.z)));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i,j+1,k),(1-w.x)*w.y*(1-w.z)));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i,j,k+1),(1-w.x)*(1-w.y)*w.z));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1,j+1,k),w.x*w.y*(1-w.z)));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1,j,k+1),w.x*(1-w.y)*w.z));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i,j+1,k+1),(1-w.x)*w.y*w.z));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1,j+1,k+1),w.x*w.y*w.z));
    return weights;
}
//#####################################################################
// Function Extrema_From_Base_Node_Helper
//#####################################################################
template<class T,class T2> VECTOR<T2,2>
Extrema_From_Base_Node_Helper(const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u_min,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u_max,const VECTOR<T,1>& X,const VECTOR<int,1>& index)
{
    int i=index.x;
    return VECTOR<T2,2>(Componentwise_Min(u_min(i),u_min(i+1)),Componentwise_Max(u_max(i),u_max(i+1)));
}
//#####################################################################
// Function Extrema_From_Base_Node_Helper
//#####################################################################
template<class T,class T2> VECTOR<T2,2>
Extrema_From_Base_Node_Helper(const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u_min,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u_max,const VECTOR<T,2>& X,const VECTOR<int,2>& index)
{
    int i=index.x,j=index.y;
    return VECTOR<T2,2>(Componentwise_Min(u_min(i,j),u_min(i+1,j),u_min(i,j+1),u_min(i+1,j+1)),Componentwise_Max(u_max(i,j),u_max(i+1,j),u_max(i,j+1),u_max(i+1,j+1)));
}
//#####################################################################
// Function Extrema_From_Base_Node_Helper
//#####################################################################
template<class T,class T2> VECTOR<T2,2>
Extrema_From_Base_Node_Helper(const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<VECTOR<T2,3> >& u_min,const ARRAYS_ND_BASE<VECTOR<T2,3> >& u_max,const VECTOR<T,3>& X,const VECTOR<int,3>& index)
{
    const T2 *min_base=&u_min(index),*max_base=&u_max(index),*min_basep1=min_base+u_min.counts.y*u_min.counts.z,*max_basep1=max_base+u_max.counts.y*u_max.counts.z;
    return VECTOR<T2,2>(
        Componentwise_Min(*min_base,*min_basep1,min_base[u_min.counts.z],min_basep1[u_min.counts.z],min_base[1],min_basep1[1],min_base[u_min.counts.z+1],min_basep1[u_min.counts.z+1]),
        Componentwise_Max(*max_base,*max_basep1,max_base[u_max.counts.z],max_basep1[u_max.counts.z],max_base[1],max_basep1[1],max_base[u_max.counts.z+1],max_basep1[u_max.counts.z+1]));
}
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& index) const
{
    return From_Base_Node_Helper(grid,u,X,index);
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> > LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node_Weights(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& index) const
{
    return From_Base_Node_Weights_Helper(grid,u,X,index);
}
//#####################################################################
// Function Extrema_From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> VECTOR<T2,2> LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Extrema_From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u_min,const T_ARRAYS_T2& u_max,const TV& X,const TV_INT& index) const
{
    return Extrema_From_Base_Node_Helper(grid,u_min,u_max,X,index);
}
//#####################################################################
// Function From_Block_Face_Component
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> typename T_GRID::VECTOR_T::SCALAR LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Block_Face_Component(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
{
    return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face_Component(axis,block,u,X);
}
//#####################################################################
// Function From_Block_Face_Component_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<FACE_INDEX<T_GRID::VECTOR_T::dimension>,typename T_GRID::VECTOR_T::SCALAR> > LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Block_Face_Component_Weights(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
{
    return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face_Component_Weights(axis,block,u,X);
}
