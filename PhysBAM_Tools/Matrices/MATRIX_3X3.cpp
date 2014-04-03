//#####################################################################
// Copyright 2002-2007, Silvia Salinas-Blemker, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Igor Neverov, Craig Schroeder, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_3X3
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Math_Tools/minabs.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> MATRIX<T,3>::
MATRIX(const MATRIX_MXN<T>& matrix_input)
{
    assert(matrix_input.n==3 && matrix_input.m==3);for(int i=0;i<3;i++) for(int j=0;j<3;j++) x[i+3*j]=(T)matrix_input(i+1,j+1);
}
//#####################################################################
// Function Higham_Iterate
//#####################################################################
template<class T> MATRIX<T,3> MATRIX<T,3>::
Higham_Iterate(const T tolerance,const int max_iterations,const bool exit_on_max_iterations) const
{
    MATRIX<T,3> X=*this;int iterations=0;
    for(;;){
        MATRIX<T,3> Y=(T).5*(X+X.Inverse_Transposed());
        if((X-Y).Max_Abs()<tolerance) return Y;
        X=Y;
        if(++iterations>=max_iterations){
            if(exit_on_max_iterations) PHYSBAM_FATAL_ERROR();
            return X;}}
}
//#####################################################################
// Function Fast_Singular_Value_Decomposition
//#####################################################################
// U and V rotations, smallest singular value possibly negative
template<class T> void MATRIX<T,3>::
Fast_Singular_Value_Decomposition(MATRIX<T,3>& U,DIAGONAL_MATRIX<T,3>& singular_values,MATRIX<T,3>& V) const // 182 mults, 112 adds, 6 divs, 11 sqrts, 1 atan2, 1 sincos
{
    if(!IS_SAME<T,double>::value){
        MATRIX<double,3> U_double,V_double;DIAGONAL_MATRIX<double,3> singular_values_double;
        MATRIX<double,3>(*this).Fast_Singular_Value_Decomposition(U_double,singular_values_double,V_double);
        U=MATRIX<T,3>(U_double);singular_values=DIAGONAL_MATRIX<T,3>(singular_values_double);V=MATRIX<T,3>(V_double);return;}
    // now T is double

    // decompose normal equations
    DIAGONAL_MATRIX<T,3> lambda;
    Normal_Equations_Matrix().Fast_Solve_Eigenproblem(lambda,V); // 18m+12a + 95m+64a+3d+5s+1atan2+1sincos

    // compute singular values
    if(lambda.x33<0) lambda=lambda.Clamp_Min(0);
    singular_values=lambda.Sqrt(); // 3s
    if(Determinant()<0) singular_values.x33=-singular_values.x33; // 9m+5a

    // compute singular vectors
    U.Column(1)=(*this*V.Column(1)).Normalized(); // 15m+8a+1d+1s
    VECTOR<T,3> v1_orthogonal=U.Column(1).Unit_Orthogonal_Vector(); // 6m+2a+1d+1s
    MATRIX<T,3,2> other_v(v1_orthogonal,VECTOR<T,3>::Cross_Product(U.Column(1),v1_orthogonal)); // 6m+3a
    U.Column(2)=other_v*(other_v.Transpose_Times(*this*V.Column(2))).Normalized(); // 6m+3a + 6m+4a + 9m+6a + 6m+2a+1d+1s = 27m+15a+1d+1s
    U.Column(3)=VECTOR<T,3>::Cross_Product(U.Column(1),U.Column(2)); // 6m+3a
}
//#####################################################################
// Function Fast_Indefinite_Polar_Decomposition
//#####################################################################
template<class T> void MATRIX<T,3>::
Fast_Indefinite_Polar_Decomposition(MATRIX<T,3>& Q,SYMMETRIC_MATRIX<T,3>& S) const
{
    MATRIX<T,3> U,V;DIAGONAL_MATRIX<T,3> D;Fast_Singular_Value_Decomposition(U,D,V);
    Q=U.Times_Transpose(V);S=SYMMETRIC_MATRIX<T,3>::Conjugate(V,D);
}
//#####################################################################
// Function Simplex_Minimum_Altitude
//#####################################################################
template<class T> T MATRIX<T,3>::
Simplex_Minimum_Altitude() const
{
    typedef VECTOR<T,3> TV;
    TV X1=Column(1),X2=Column(2),X3=Column(3);
    return minabs(
        TV::Dot_Product(X1,TV::Cross_Product(X2-X1,X3-X1).Normalized()),
        TV::Dot_Product(X2-X1,TV::Cross_Product(X3,X2).Normalized()),
        TV::Dot_Product(X3-X2,TV::Cross_Product(X1,X3).Normalized()),
        TV::Dot_Product(X3,TV::Cross_Product(X1,X2).Normalized()));
}
//#####################################################################
// Function Componentwise_Min
//#####################################################################
template<class T> MATRIX<T,3> MATRIX<T,3>::
Componentwise_Min(const MATRIX& v1,const MATRIX& v2)
{
    return MATRIX(min(v1.x[0],v2.x[0]),min(v1.x[1],v2.x[1]),min(v1.x[2],v2.x[2]),min(v1.x[3],v2.x[3]),min(v1.x[4],v2.x[4]),min(v1.x[5],v2.x[5]),
        min(v1.x[6],v2.x[6]),min(v1.x[7],v2.x[7]),min(v1.x[8],v2.x[8]));
}
//#####################################################################
// Function Componentwise_Max
//#####################################################################
template<class T> MATRIX<T,3> MATRIX<T,3>::
Componentwise_Max(const MATRIX& v1,const MATRIX& v2)
{
    return MATRIX(max(v1.x[0],v2.x[0]),max(v1.x[1],v2.x[1]),max(v1.x[2],v2.x[2]),max(v1.x[3],v2.x[3]),max(v1.x[4],v2.x[4]),max(v1.x[5],v2.x[5]),
        max(v1.x[6],v2.x[6]),max(v1.x[7],v2.x[7]),max(v1.x[8],v2.x[8]));
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> MATRIX_MXN<T> MATRIX<T,3>::
operator*(const MATRIX_MXN<T>& A) const
{
    assert(3==A.m);
    MATRIX_MXN<T> matrix(3,A.n);
    for(int j=1;j<=A.n;j++) for(int i=1;i<=3;i++) for(int k=1;k<=3;k++) matrix(i,j)+=(*this)(i,k)*A(k,j);
    return matrix;
}
//#####################################################################
// Function Inverse
//#####################################################################
template<class T> MATRIX<T,3> MATRIX<T,3>::
Inverse() const
{
    T cofactor11=x[4]*x[8]-x[7]*x[5],cofactor12=x[7]*x[2]-x[1]*x[8],cofactor13=x[1]*x[5]-x[4]*x[2];
    T determinant=x[0]*cofactor11+x[3]*cofactor12+x[6]*cofactor13;
    assert(determinant!=0);
    T s=1/determinant;
    return s*MATRIX(cofactor11,cofactor12,cofactor13,x[6]*x[5]-x[3]*x[8],x[0]*x[8]-x[6]*x[2],x[3]*x[2]-x[0]*x[5],x[3]*x[7]-x[6]*x[4],x[6]*x[1]-x[0]*x[7],x[0]*x[4]-x[3]*x[1]);
}
//#####################################################################
// Function Solve_Linear_System
//#####################################################################
template<class T> VECTOR<T,3> MATRIX<T,3>::
Solve_Linear_System(const VECTOR<T,3>& b) const // 33 mults, 17 adds, 1 div
{
    T cofactor11=x[4]*x[8]-x[7]*x[5],cofactor12=x[7]*x[2]-x[1]*x[8],cofactor13=x[1]*x[5]-x[4]*x[2];
    T determinant=x[0]*cofactor11+x[3]*cofactor12+x[6]*cofactor13;
    assert(determinant!=0);
    return MATRIX(cofactor11,cofactor12,cofactor13,x[6]*x[5]-x[3]*x[8],x[0]*x[8]-x[6]*x[2],x[3]*x[2]-x[0]*x[5],x[3]*x[7]-x[6]*x[4],x[6]*x[1]-x[0]*x[7],x[0]*x[4]-x[3]*x[1])*b/determinant;
}
//#####################################################################
// Function Robust_Solve_Linear_System
//#####################################################################
template<class T> VECTOR<T,3> MATRIX<T,3>::
Robust_Solve_Linear_System(const VECTOR<T,3>& b) const // 34 mults, 17 adds, 1 div
{
    T cofactor11=x[4]*x[8]-x[7]*x[5],cofactor12=x[7]*x[2]-x[1]*x[8],cofactor13=x[1]*x[5]-x[4]*x[2];
    T determinant=x[0]*cofactor11+x[3]*cofactor12+x[6]*cofactor13;
    VECTOR<T,3> unscaled_result=MATRIX(cofactor11,cofactor12,cofactor13,x[6]*x[5]-x[3]*x[8],x[0]*x[8]-x[6]*x[2],x[3]*x[2]-x[0]*x[5],x[3]*x[7]-x[6]*x[4],x[6]*x[1]-x[0]*x[7],
        x[0]*x[4]-x[3]*x[1])*b;
    T relative_tolerance=(T)FLT_MIN*unscaled_result.Max_Abs();
    if(abs(determinant)<=relative_tolerance){relative_tolerance=max(relative_tolerance,(T)FLT_MIN);determinant=determinant>=0?relative_tolerance:-relative_tolerance;}
    return unscaled_result/determinant;
}
//#####################################################################
// Function Q_From_QR_Factorization
//#####################################################################
template<class T> MATRIX<T,3> MATRIX<T,3>::
Q_From_QR_Factorization() const // Gram Schmidt
{
    int k;MATRIX Q=*this;
    T one_over_r11=1/sqrt((sqr(Q.x[0])+sqr(Q.x[1])+sqr(Q.x[2])));
    for(k=0;k<=2;k++) Q.x[k]=one_over_r11*Q.x[k];
    T r12=Q.x[0]*Q.x[3]+Q.x[1]*Q.x[4]+Q.x[2]*Q.x[5];
    Q.x[3]-=r12*Q.x[0];
    Q.x[4]-=r12*Q.x[1];
    Q.x[5]-=r12*Q.x[2];
    T r13=Q.x[0]*Q.x[6]+Q.x[1]*Q.x[7]+Q.x[2]*Q.x[8];
    Q.x[6]-=r13*Q.x[0];
    Q.x[7]-=r13*Q.x[1];
    Q.x[8]-=r13*Q.x[2];
    T one_over_r22=1/sqrt((sqr(Q.x[3])+sqr(Q.x[4])+sqr(Q.x[5])));
    for(k=3;k<=5;k++) Q.x[k]=one_over_r22*Q.x[k];
    T r23=Q.x[3]*Q.x[6]+Q.x[4]*Q.x[7]+Q.x[5]*Q.x[8];
    Q.x[6]-=r23*Q.x[3];
    Q.x[7]-=r23*Q.x[4];
    Q.x[8]-=r23*Q.x[5];
    T one_over_r33=1/sqrt((sqr(Q.x[6])+sqr(Q.x[7])+sqr(Q.x[8])));
    for(k=6;k<=8;k++) Q.x[k]=one_over_r33*Q.x[k];
    return Q;
}
//#####################################################################
// Function R_From_QR_Factorization
//#####################################################################
template<class T> UPPER_TRIANGULAR_MATRIX<T,3> MATRIX<T,3>::
R_From_QR_Factorization() const // Gram Schmidt
{
    int k;
    MATRIX Q=*this;
    UPPER_TRIANGULAR_MATRIX<T,3> R;
    R.x11=sqrt((sqr(Q.x[0])+sqr(Q.x[1])+sqr(Q.x[2])));
    T one_over_r11=1/R.x11;
    for(k=0;k<=2;k++) Q.x[k]=one_over_r11*Q.x[k];
    R.x12=Q.x[0]*Q.x[3]+Q.x[1]*Q.x[4]+Q.x[2]*Q.x[5];
    Q.x[3]-=R.x12*Q.x[0];
    Q.x[4]-=R.x12*Q.x[1];
    Q.x[5]-=R.x12*Q.x[2];
    R.x13=Q.x[0]*Q.x[6]+Q.x[1]*Q.x[7]+Q.x[2]*Q.x[8];
    Q.x[6]-=R.x13*Q.x[0];
    Q.x[7]-=R.x13*Q.x[1];
    Q.x[8]-=R.x13*Q.x[2];
    R.x22=sqrt((sqr(Q.x[3])+sqr(Q.x[4])+sqr(Q.x[5])));
    T one_over_r22=1/R.x22;
    for(k=3;k<=5;k++) Q.x[k]=one_over_r22*Q.x[k];
    R.x23=Q.x[3]*Q.x[6]+Q.x[4]*Q.x[7]+Q.x[5]*Q.x[8];
    Q.x[6]-=R.x23*Q.x[3];
    Q.x[7]-=R.x23*Q.x[4];
    Q.x[8]-=R.x23*Q.x[5];
    R.x33=sqrt((sqr(Q.x[6])+sqr(Q.x[7])+sqr(Q.x[8])));
    return R;
}
//#####################################################################
template class MATRIX<float,3>;
template class MATRIX<double,3>;
