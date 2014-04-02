//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BANDED_SYMMETRIC_MATRIX
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Math_Tools/minabs.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Matrices/BANDED_SYMMETRIC_MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <limits>
namespace PhysBAM{
//#####################################################################
// Function Multiply
//#####################################################################
template<class T,int bandwidth> void BANDED_SYMMETRIC_MATRIX<T,bandwidth>::
Multiply(const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const
{
    STATIC_ASSERT(bandwidth==3); // TODO: implement for non-tridiagonal matrices
    assert(A.Size()==Size() && x.Size()==Size());
    T x_previous=x(1);
    result(1)=A(1)[1]*x_previous;
    for(int i=2;i<=Size();i++){
        T x_i=x(i);
        result(i-1)+=A(i-1)[2]*x_i;
        result(i)=A(i-1)[2]*x_previous+A(i)[1]*x_i;
        x_previous=x_i;}
}
//#####################################################################
// Function Power_Iterate_Shifted
//#####################################################################
template<class T,int bandwidth> bool BANDED_SYMMETRIC_MATRIX<T,bandwidth>::
Power_Iterate_Shifted(VECTOR_ND<T>& x,const T shift,T& eigenvalue,const T tolerance,const int max_iterations) const
{
    VECTOR_ND<T> y(Size(),false);
    T tolerance_squared_over_2=sqr(tolerance)/2;
    for(int iteration=1;iteration<=max_iterations;iteration++){
        T magnitude=x.Magnitude();if(!magnitude) return false;
        x/=magnitude;
        Multiply(x,y);if(shift) y+=shift*x;
        eigenvalue=VECTOR_ND<T>::Dot_Product(x,y);
        x=y-eigenvalue*x;
        T error_squared=x.Magnitude_Squared();
        x=y;
        if(error_squared<=tolerance_squared_over_2) return true;}
    return false;
}
//#####################################################################
// Function Eigenvalue_Range
//#####################################################################
template<class T,int bandwidth> bool BANDED_SYMMETRIC_MATRIX<T,bandwidth>::
Eigenvalue_Range(VECTOR_ND<T>& x,INTERVAL<T>& eigenvalue_range,const T tolerance,const int max_iterations) const
{
    if(Size()==1){eigenvalue_range=INTERVAL<T>(A(1)(1));return true;}
    T lambda_min,lambda_max;
    if(!Power_Iterate(x,lambda_max,tolerance,max_iterations)) return false;
    eigenvalue_range=INTERVAL<T>(lambda_max);
    x.Set_To_Orthogonal_Vector(); // change vector to hopefully include component along other eigenvalue
    if(!Power_Iterate_Shifted(x,-lambda_max,lambda_min,tolerance,max_iterations)) return false;
    eigenvalue_range.Enlarge_Nonempty_Box_To_Include_Point(lambda_min+lambda_max);
    return true;
}
//#####################################################################
// Function Implicit_QR_Step
//#####################################################################
template<class T,int bandwidth> void BANDED_SYMMETRIC_MATRIX<T,bandwidth>::
Implicit_QR_Step(const int start,const T shift)
{
    // see Golub and Van Loan, 8.3.2, p. 420
    T x=A(start)[1]-shift,z=A(start)[2];
    for(int k=start;;k++){ // chase offdiagonal element z into oblivion
        T c,s;Givens(x,z,c,s);
        if(k>start){T& b=A(k-1)[2];b=c*b-s*z;}
        T &ap=A(k)[1],&bp=A(k)[2],&aq=A(k+1)[1],&bq=A(k+1)[2]; // p,q = k,k+1
        T cc=c*c,cs=c*s,ss=s*s,ccap=cc*ap,csbp=cs*bp,ssaq=ss*aq,csap=cs*ap,ccbp=cc*bp,ssbp=ss*bp,csaq=cs*aq,ssap=ss*ap,ccaq=cc*aq;
        ap=ccap-2*csbp+ssaq;bp=csap+ccbp-ssbp-csaq;aq=ssap+2*csbp+ccaq; // TODO: save 7 flops somehow
        if(!bq) return; // stop if matrix decouples
        x=bp;z=-s*bq;bq=c*bq;}
}
//#####################################################################
// Function Diagonalize
//#####################################################################
template<class T,int bandwidth> void BANDED_SYMMETRIC_MATRIX<T,bandwidth>::
Diagonalize()
{
    static const T tolerance=10*std::numeric_limits<T>::epsilon();
    int m=Size();if(m<2) return;
    A(m)[2]=0; // use a_m,m+1 as a sentinel
    int start=1;
    for(;;){
        // find unreduced section
        while(!A(start)[2] || abs(A(start)[2])<=tolerance*maxabs(A(start)[1],A(start)[2])){
            A(start)[2]=0;if(++start==m) return;}
        int end=start+1;
        while(A(end)[2] && abs(A(end)[2])>tolerance*maxabs(A(end)[1],A(end)[2])) end++;
        A(end)[2]=0;
        // compute shift
        T ap=A(end-1)[1],bp=A(end-1)[2],aq=A(end)[1],d=(ap-aq)/2;
        T shift=d?aq-bp*bp/(d+sign(d)*sqrt(d*d+bp*bp)):aq-bp;
        Implicit_QR_Step(start,shift);}
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
template<class T,int bandwidth> void BANDED_SYMMETRIC_MATRIX<T,bandwidth>::
Eigenvalues(ARRAY<T>& eigenvalues) const
{
    BANDED_SYMMETRIC_MATRIX D(*this);
    D.Diagonalize();
    eigenvalues.Resize(Size(),false,false);
    for(int i=1;i<=Size();i++) eigenvalues(i)=D.A(i)[1];
}
//#####################################################################
// Function Print_Spectral_Information
//#####################################################################
template<class T,int bandwidth> void BANDED_SYMMETRIC_MATRIX<T,bandwidth>::
Print_Spectral_Information() const
{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(!Size()){
        LOG::cout<<"eigenvalue range = null, condition = null"<<std::endl;
        return;}
    ARRAY<T> D;Eigenvalues(D);
    T lambda_min=ARRAYS_COMPUTATIONS::Min(D),lambda_max=ARRAYS_COMPUTATIONS::Max(D);
    T condition=lambda_min*lambda_max>0?Robust_Divide(maxabs(lambda_min,lambda_max),minabs(lambda_min,lambda_max)):0;
    LOG::cout<<"eigenvalue range = "<<lambda_min<<" "<<lambda_max<<", condition = "<<condition<<std::endl;
    Sort(D);LOG::cout<<"eigenvalues = "<<D;
#endif
}
//#####################################################################
template class BANDED_SYMMETRIC_MATRIX<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BANDED_SYMMETRIC_MATRIX<double,3>;
#endif
}
