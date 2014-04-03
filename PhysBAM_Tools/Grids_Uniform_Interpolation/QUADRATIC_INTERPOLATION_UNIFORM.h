//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUADRATIC_INTERPOLATION_UNIFORM 
//#####################################################################
#ifndef __QUADRATIC_INTERPOLATION_UNIFORM__
#define __QUADRATIC_INTERPOLATION_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<T_GRID>
class QUADRATIC_INTERPOLATION_UNIFORM:public INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    //typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> > T_ARRAYS_T2;
public:
    typedef INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP> BASE;
    using BASE::Clamped_Index_Interior_End_Minus_One;

    int a_scheme;
public:

    QUADRATIC_INTERPOLATION_UNIFORM();
    ~QUADRATIC_INTERPOLATION_UNIFORM();

    T2 Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE;
    ARRAY<PAIR<TV_INT,T> > Clamped_To_Array_Weights(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE;
    T2 From_Base_Node_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const;
    //Sanity Check
    T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const;

    //New method averaging
    /*T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const
    {T w=(X.x-grid.X(index.x).x)*grid.one_over_dX.x;T a=0;
    if(a_scheme==1) a=(u(index.x-1)+u(index.x+2)-u(index.x)-u(index.x+1))/4.;
    else if(a_scheme==2){T a1=(u(index.x+2)+u(index.x)-2*u(index.x+1))/2.;T a2=(u(index.x+1)+u(index.x-1)-2*u(index.x))/2.;a=minmag(a1,a2);}
    return a*w*(w-1)+(1-w)*u(index.x)+w*u(index.x+1);}*/

    template<class T,int d>
    T abs(const VECTOR<T,d>& vector) const
    {return vector.Magnitude();}
    
    template<class T>
    T abs(const T val) const
    {return val>=0?val:-1*val;}

    //New method weights
    ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const;

    //Weights with alpha
    /*ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const
    {ARRAY<PAIR<TV_INT,T> > weights;T w=(X.x-grid.X(index.x).x)*grid.one_over_dX.x;
    T2 phi=From_Base_Node_Helper(grid,u,X,index);T2 phi_linear=(1-w)*u(index.x)+w*u(index.x+1);T alpha=(u(index.x)+u(index.x+1))?(phi-phi_linear)/(u(index.x)+u(index.x+1)):0;
    weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),1-w+alpha));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w+alpha));
    for(int i=1;i<=weights.m;i++) weights(i).y=max((T)0.,min((T)1.,weights(i).y));
    T error=phi-(1-w+alpha)*u(index.x)-(w+alpha)*u(index.x+1);
    if(abs(error)>1e-5) std::cout<<"Error is "<<error<<" and weights are "<<weights<<std::endl;
    return weights;}*/

    //Linear Programming
    /*ARRAY<PAIR<TV_INT,T> > From_Base_Node_W GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const
    {ARRAY<PAIR<TV_INT,T> > weights;T w=(X.x-grid.X(index.x).x)*grid.one_over_dX.x;
    T2 phi=From_Base_Node_Helper(grid,u,X,index);T2 phi_linear=(1-w)*u(index.x)+w*u(index.x+1);T alpha=(u(index.x)+u(index.x+1))?(phi-phi_linear)/(u(index.x)+u(index.x+1)):0;
    weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),1-w+alpha));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w+alpha));
    for(int i=1;i<=weights.m;i++) weights(i).y=max((T)0.,min((T)1.,weights(i).y));
    T error=phi-(1-w+alpha)*u(index.x)-(w+alpha)*u(index.x+1);
    if(abs(error)>1e-5) std::cout<<"Error is "<eights(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const
    {ARRAY<PAIR<TV_INT,T> > weights;T2 phi=From_Base_Node(grid,u,X,index);T w=(X.x-grid.X(index.x).x)*grid.one_over_dX.x;
    T leftDxx=(u(index.x+1)-2*u(index.x)+u(index.x-1)),rightDxx=(u(index.x+2)-2*u(index.x+1)+u(index.x));
    MATRIX_MXN<T> A(2,4);A(1,1)=1;A(2,1)=0;A(2,2)=1;A(2,3)=1;A(2,4)=1;VECTOR_ND<T> x(4);VECTOR_ND<T> b(2);b(2)=1;
    if(abs(leftDxx)<abs(rightDxx)){
        A(1,2)=-u(index.x-1);A(1,3)=-u(index.x);A(1,4)=-u(index.x+1);b(1)=-phi;
        x(2)=0;x(3)=1-w;x(4)=w;x(1)=u(index.x-1)*x(2)+u(index.x)*x(3)+u(index.x+1)*x(4);
        weights.Append(PAIR<TV_INT,T>(TV_INT(index.x-1),0));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),0));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),0));}
    else{
        A(1,2)=-u(index.x);A(1,3)=-u(index.x+1);A(1,4)=-u(index.x+2);b(1)=-phi;
        x(2)=1-w;x(3)=w;x(4)=0;x(1)=u(index.x)*x(2)+u(index.x+1)*x(3)+u(index.x+2)*x(4);
        weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),0));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),0));weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+2),0));}
    //Calculate x
    for(int i=1;i<=3;i++) weights(i).y=x(i+1);
    return weights;}*/
    
    //Clamp Version of Gibou
    /*ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const
    {ARRAY<PAIR<TV_INT,T> > weights;T w=(X.x-grid.X(index.x).x)*grid.one_over_dX.x;T one_over_dX2=grid.one_over_dX.x*grid.one_over_dX.x;
    T leftDxx=(u(index.x+1)-2*u(index.x)+u(index.x-1)),rightDxx=(u(index.x+2)-2*u(index.x+1)+u(index.x));
    if(abs(leftDxx)<abs(rightDxx)){
        weights.Append(PAIR<TV_INT,T>(TV_INT(index.x-1),-one_over_dX2*(w*(1-w))/2.));
        weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),(1-w)+2*one_over_dX2*(w*(1-w))/2.));
        weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w-one_over_dX2*(w*(1-w))/2.));}
    else{
        weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),(1-w)-one_over_dX2*(w*(1-w))/2.));
        weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w+2*one_over_dX2*(w*(1-w))/2.));
        weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+2),-one_over_dX2*(w*(1-w))/2.));}
    //weights.Append(PAIR<TV_INT,T>(TV_INT(index.x-1),w*(w-1)/2.));
    //weights.Append(PAIR<TV_INT,T>(TV_INT(index.x),w*w-1));
    //weights.Append(PAIR<TV_INT,T>(TV_INT(index.x+1),w*(w+1)/2.));
    for(int i=1;i<=weights.m;i++) weights(i).y=max((T)0.,min((T)1.,weights(i).y));
    return weights;}*/

    T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const;
    ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const;
    T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const;
    ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const;

//#####################################################################
};
}
#endif
