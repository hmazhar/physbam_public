//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_4X4
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
using namespace PhysBAM;
//#####################################################################
// Function operator*
//#####################################################################
template<class T> MATRIX<T,4> MATRIX<T,4>::
operator*(const MATRIX& A) const
{
    return MATRIX(x[0]*A.x[0]+x[4]*A.x[1]+x[8]*A.x[2]+x[12]*A.x[3],x[1]*A.x[0]+x[5]*A.x[1]+x[9]*A.x[2]+x[13]*A.x[3],
        x[2]*A.x[0]+x[6]*A.x[1]+x[10]*A.x[2]+x[14]*A.x[3],x[3]*A.x[0]+x[7]*A.x[1]+x[11]*A.x[2]+x[15]*A.x[3],
        x[0]*A.x[4]+x[4]*A.x[5]+x[8]*A.x[6]+x[12]*A.x[7],x[1]*A.x[4]+x[5]*A.x[5]+x[9]*A.x[6]+x[13]*A.x[7],
        x[2]*A.x[4]+x[6]*A.x[5]+x[10]*A.x[6]+x[14]*A.x[7],x[3]*A.x[4]+x[7]*A.x[5]+x[11]*A.x[6]+x[15]*A.x[7],
        x[0]*A.x[8]+x[4]*A.x[9]+x[8]*A.x[10]+x[12]*A.x[11],x[1]*A.x[8]+x[5]*A.x[9]+x[9]*A.x[10]+x[13]*A.x[11],
        x[2]*A.x[8]+x[6]*A.x[9]+x[10]*A.x[10]+x[14]*A.x[11],x[3]*A.x[8]+x[7]*A.x[9]+x[11]*A.x[10]+x[15]*A.x[11],
        x[0]*A.x[12]+x[4]*A.x[13]+x[8]*A.x[14]+x[12]*A.x[15],x[1]*A.x[12]+x[5]*A.x[13]+x[9]*A.x[14]+x[13]*A.x[15],
        x[2]*A.x[12]+x[6]*A.x[13]+x[10]*A.x[14]+x[14]*A.x[15],x[3]*A.x[12]+x[7]*A.x[13]+x[11]*A.x[14]+x[15]*A.x[15]);
}
//#####################################################################
// Function Inverse
//#####################################################################
template<class T> MATRIX<T,4> MATRIX<T,4>::
Inverse() const
{
    int pivot_row;
    T p1=abs(x[0]),p2=abs(x[1]),p3=abs(x[2]),p4=abs(x[3]);
    if(p1>p2){if(p1>p3){if(p1>p4) pivot_row=1;else pivot_row=4;}else if(p3>p4) pivot_row=3;else pivot_row=4;}
    else if(p2>p3){if(p2>p4) pivot_row=2;else pivot_row=4;}
    else if(p3>p4) pivot_row=3;else pivot_row=4;
    T a=x[pivot_row-1];assert(a!=0);
    VECTOR<T,3> b(x[3+pivot_row],x[7+pivot_row],x[11+pivot_row]),c;MATRIX<T,3> d;
    switch(pivot_row){
        case 1:c=VECTOR<T,3>(x[1],x[2],x[3]);d=MATRIX<T,3>(x[5],x[6],x[7],x[9],x[10],x[11],x[13],x[14],x[15]);break;
        case 2:c=VECTOR<T,3>(x[0],x[2],x[3]);d=MATRIX<T,3>(x[4],x[6],x[7],x[8],x[10],x[11],x[12],x[14],x[15]);break;
        case 3:c=VECTOR<T,3>(x[0],x[1],x[3]);d=MATRIX<T,3>(x[4],x[5],x[7],x[8],x[9],x[11],x[12],x[13],x[15]);break;
        default:c=VECTOR<T,3>(x[0],x[1],x[2]);d=MATRIX<T,3>(x[4],x[5],x[6],x[8],x[9],x[10],x[12],x[13],x[14]);}
    T m=-1/a;b*=m;d+=MATRIX<T,3>::Outer_Product(c,b); // find schur complement
    MATRIX<T,3> h=d.Inverse();VECTOR<T,3> g=h*(c*m);VECTOR<T,3> f=b*h;T e=VECTOR<T,3>::Dot_Product(b,g)-m; // compute parts of inverse
    switch(pivot_row){
        case 1:return MATRIX(e,g.x,g.y,g.z,f.x,h(1,1),h(2,1),h(3,1),f.y,h(1,2),h(2,2),h(3,2),f.z,h(1,3),h(2,3),h(3,3));break;
        case 2:return MATRIX(f.x,h(1,1),h(2,1),h(3,1),e,g.x,g.y,g.z,f.y,h(1,2),h(2,2),h(3,2),f.z,h(1,3),h(2,3),h(3,3));break;
        case 3:return MATRIX(f.x,h(1,1),h(2,1),h(3,1),f.y,h(1,2),h(2,2),h(3,2),e,g.x,g.y,g.z,f.z,h(1,3),h(2,3),h(3,3));break;}
    return MATRIX(f.x,h(1,1),h(2,1),h(3,1),f.y,h(1,2),h(2,2),h(3,2),f.z,h(1,3),h(2,3),h(3,3),e,g.x,g.y,g.z);
}
//#####################################################################
// Function Cofactor_Matrix
//#####################################################################
template<class T> MATRIX<T,4> MATRIX<T,4>::
Cofactor_Matrix() const
{
    T y[16];
    y[0]=x[5]*x[10]*x[15]+x[9]*x[14]*x[7]+x[13]*x[6]*x[11]-x[5]*x[14]*x[11]-x[9]*x[6]*x[15]-x[13]*x[10]*x[7];
    y[1]=-x[6]*x[11]*x[12]-x[10]*x[15]*x[4]-x[14]*x[7]*x[8]+x[6]*x[15]*x[8]+x[10]*x[7]*x[12]+x[14]*x[11]*x[4];
    y[2]=x[7]*x[8]*x[13]+x[11]*x[12]*x[5]+x[15]*x[4]*x[9]-x[7]*x[12]*x[9]-x[11]*x[4]*x[13]-x[15]*x[8]*x[5];
    y[3]=-x[4]*x[9]*x[14]-x[8]*x[13]*x[6]-x[12]*x[5]*x[10]+x[4]*x[13]*x[10]+x[8]*x[5]*x[14]+x[12]*x[9]*x[6];
    y[4]=-x[9]*x[14]*x[3]-x[13]*x[2]*x[11]-x[1]*x[10]*x[15]+x[9]*x[2]*x[15]+x[13]*x[10]*x[3]+x[1]*x[14]*x[11];
    y[5]=x[10]*x[15]*x[0]+x[14]*x[3]*x[8]+x[2]*x[11]*x[12]-x[10]*x[3]*x[12]-x[14]*x[11]*x[0]-x[2]*x[15]*x[8];
    y[6]=-x[11]*x[12]*x[1]-x[15]*x[0]*x[9]-x[3]*x[8]*x[13]+x[11]*x[0]*x[13]+x[15]*x[8]*x[1]+x[3]*x[12]*x[9];
    y[7]=x[8]*x[13]*x[2]+x[12]*x[1]*x[10]+x[0]*x[9]*x[14]-x[8]*x[1]*x[14]-x[12]*x[9]*x[2]-x[0]*x[13]*x[10];
    y[8]=x[13]*x[2]*x[7]+x[1]*x[6]*x[15]+x[5]*x[14]*x[3]-x[13]*x[6]*x[3]-x[1]*x[14]*x[7]-x[5]*x[2]*x[15];
    y[9]=-x[14]*x[3]*x[4]-x[2]*x[7]*x[12]-x[6]*x[15]*x[0]+x[14]*x[7]*x[0]+x[2]*x[15]*x[4]+x[6]*x[3]*x[12];
    y[10]=x[15]*x[0]*x[5]+x[3]*x[4]*x[13]+x[7]*x[12]*x[1]-x[15]*x[4]*x[1]-x[3]*x[12]*x[5]-x[7]*x[0]*x[13];
    y[11]=-x[12]*x[1]*x[6]-x[0]*x[5]*x[14]-x[4]*x[13]*x[2]+x[12]*x[5]*x[2]+x[0]*x[13]*x[6]+x[4]*x[1]*x[14];
    y[12]=-x[1]*x[6]*x[11]-x[5]*x[10]*x[3]-x[9]*x[2]*x[7]+x[1]*x[10]*x[7]+x[5]*x[2]*x[11]+x[9]*x[6]*x[3];
    y[13]=x[2]*x[7]*x[8]+x[6]*x[11]*x[0]+x[10]*x[3]*x[4]-x[2]*x[11]*x[4]-x[6]*x[3]*x[8]-x[10]*x[7]*x[0];
    y[14]=-x[3]*x[4]*x[9]-x[7]*x[8]*x[1]-x[11]*x[0]*x[5]+x[3]*x[8]*x[5]+x[7]*x[0]*x[9]+x[11]*x[4]*x[1];
    y[15]=x[0]*x[5]*x[10]+x[4]*x[9]*x[2]+x[8]*x[1]*x[6]-x[0]*x[9]*x[6]-x[4]*x[1]*x[10]-x[8]*x[5]*x[2];
    return MATRIX(y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10],y[11],y[12],y[13],y[14],y[15]);
}
//#####################################################################
// Function Frobenius_Norm_Squared
//#####################################################################
template<class T> T MATRIX<T,4>::
Frobenius_Norm_Squared() const
{
    T sum=0;
    for(int i=0;i<16;i++) sum+=sqr(x[i]);
    return sum;
}
//#####################################################################
template class MATRIX<float,4>;
template class MATRIX<double,4>;
