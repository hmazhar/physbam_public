//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
QUADRATIC_INTERPOLATION_UNIFORM()
    :a_scheme(1)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
~QUADRATIC_INTERPOLATION_UNIFORM()
{
}
//#####################################################################
// Function Clamped_To_Array
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE
{
    return From_Base_Node(grid,u,X,INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::Clamped_Index_Interior_End_Minus_One(grid,u,X));
}
//#####################################################################
// Function Clamped_To_Array_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> > QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Clamped_To_Array_Weights(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE
{
    return From_Base_Node_Weights(grid,u,X,INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::Clamped_Index_Interior_End_Minus_One(grid,u,X));
}
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const
{
    T w=(X.x-grid.X(index.x).x)*grid.one_over_dX.x;T one_over_dX2=1;//grid.one_over_dX.x*grid.one_over_dX.x;
    T leftDxx=(u(index.x+1)-2*u(index.x)+u(index.x-1)),rightDxx=(u(index.x+2)-2*u(index.x+1)+u(index.x));
    if(abs(leftDxx)<abs(rightDxx)) return u(index.x)*(1-w)+u(index.x+1)*w-leftDxx*one_over_dX2*w*(1-w)/2.;
    else return u(index.x)*(1-w)+u(index.x+1)*w-rightDxx*one_over_dX2*w*(1-w)/2.;
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const
{
    ARRAY<PAIR<TV_INT,T> > weights=From_Base_Node_Weights(grid,u,X,index);
    T2 sum=T2();
    for(int i=1;i<=weights.m;i++) sum+=u(weights(i).x)*weights(i).y;
    return sum;
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> > QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const
{
    ARRAY<PAIR<TV_INT,T> > weights;
    T w=(X.x-grid.X(index.x).x)*grid.one_over_dX.x;
    if(a_scheme==1){
        if(abs(u(index.x))>1e-5 && abs(u(index.x+1))>1e-5){weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),(T)1-w-w*(w-(T)1)/(T)4*((T)1-u(index.x-1)/u(index.x))));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w-w*(w-(T)1)/(T)4*((T)1-u(index.x+2)/u(index.x+1))));}
        else if(abs(u(index.x))>1e-5){weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),(T)1-w-w*(w-(T)1)/(T)4*((T)1-(u(index.x-1)-u(index.x+2))/u(index.x))));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w-w*(w-(T)1)/(T)4));}
        else if(abs(u(index.x+1))>1e-5){weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),(T)1-w-w*(w-(T)1)/(T)4));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w-w*(w-(T)1)/(T)4*((T)1-(u(index.x-1)-u(index.x+2))/u(index.x+1))));}
        else{weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),(T)1-w));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w));}}
    else if(a_scheme==2){
        T a1=(u(index.x+2)+u(index.x)-(T)2*u(index.x+1))/(T)2;T a2=(u(index.x+1)+u(index.x-1)-(T)2*u(index.x))/(T)2;
        if(abs(a1)<abs(a2)){
            if(abs(u(index.x+1))>1e-5){weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),(T)1-w+w*(w-(T)1)/(T)2));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w-w*(w-(T)1)*((T)1-u(index.x+2)/((T)2*u(index.x+1)))));}
            else if(abs(u(index.x))>1e-5){weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),(T)1-w+w*(w-(T)1)/(T)2*(1.+u(index.x+2)/u(index.x))));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w-w*(w-(T)1)));}
            else{weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),(T)1-w));weights.Append(PAIR<TV_INT,T2>(TV_INT(index.x+1),w));}}
        else{
            if(abs(u(index.x))>1e-5){weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),(T)1-w-w*(w-(T)1)*((T)1-u(index.x-1)/((T)2*u(index.x)))));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w+w*(w-(T)1)/(T)2));}
            else if(abs(u(index.x+1))>1e-5){weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),(T)1-w-w*(w-(T)1)));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w+w*(w-(T)1)/(T)2*((T)1+u(index.x-1)/(u(index.x+1)))));}
            else{weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),(T)1-w));weights.Append(PAIR<TV_INT,T2>(TV_INT(index.x+1),w));}}}
    return weights;
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const
{
    PHYSBAM_NOT_IMPLEMENTED();
    return T2();
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> > QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const
{
    PHYSBAM_NOT_IMPLEMENTED();
    return ARRAY<PAIR<TV_INT,T> >();
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const
{
    PHYSBAM_NOT_IMPLEMENTED();
    return T2();
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> > QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const
{
    PHYSBAM_NOT_IMPLEMENTED();
    return ARRAY<PAIR<TV_INT,T> >();
}

template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >::QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >::~QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::~QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::~QUADRATIC_INTERPOLATION_UNIFORM();
template float QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >::Clamped_To_Array(GRID<VECTOR<float,1> > const&,ARRAYS_ND_BASE<VECTOR<float,1> > const&,VECTOR<float,1> const&) const;
template float QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::Clamped_To_Array(GRID<VECTOR<float,2> > const&,ARRAYS_ND_BASE<VECTOR<float,2> > const&,VECTOR<float,2> const&) const;
template float QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::Clamped_To_Array(GRID<VECTOR<float,3> > const&,ARRAYS_ND_BASE<VECTOR<float,3> > const&,VECTOR<float,3> const&) const;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >::QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::~QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::~QUADRATIC_INTERPOLATION_UNIFORM();
template double QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::Clamped_To_Array(GRID<VECTOR<double,2> > const&,ARRAYS_ND_BASE<VECTOR<double,2> > const&,VECTOR<double,2> const&) const;
template double QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::Clamped_To_Array(GRID<VECTOR<double,3> > const&,ARRAYS_ND_BASE<VECTOR<double,3> > const&,VECTOR<double,3> const&) const;
#endif
