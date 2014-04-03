//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_3X2
//#####################################################################
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
using namespace PhysBAM;
//#####################################################################
// Function Fast_Singular_Value_Decomposition
//#####################################################################
template<class T> void MATRIX<T,3,2>::
Fast_Singular_Value_Decomposition(MATRIX<T,3,2>& U,DIAGONAL_MATRIX<T,2>& singular_values,MATRIX<T,2>& V) const
{
    if(!IS_SAME<T,double>::value){
        MATRIX<double,3,2> U_double;DIAGONAL_MATRIX<double,2> singular_values_double;MATRIX<double,2> V_double;
        MATRIX<double,3,2>(*this).Fast_Singular_Value_Decomposition(U_double,singular_values_double,V_double);
        U=MATRIX<T,3,2>(U_double);singular_values=DIAGONAL_MATRIX<T,2>(singular_values_double);V=MATRIX<T,2>(V_double);return;}
    // now T is double

    DIAGONAL_MATRIX<T,2> lambda;Normal_Equations_Matrix().Solve_Eigenproblem(lambda,V);
    if(lambda.x22<0) lambda=lambda.Clamp_Min(0);
    singular_values=lambda.Sqrt();
    U.Column(1)=(*this*V.Column(1)).Normalized();
    VECTOR<T,3> other=VECTOR<T,3>::Cross_Product(Weighted_Normal(),U.Column(1));
    T other_magnitude=other.Magnitude();
    U.Column(2)=other_magnitude?other/other_magnitude:U.Column(1).Unit_Orthogonal_Vector();
}
//#####################################################################
// Function Fast_Indefinite_Polar_Decomposition
//#####################################################################
template<class T> void MATRIX<T,3,2>::
Fast_Indefinite_Polar_Decomposition(MATRIX<T,3,2>& Q,SYMMETRIC_MATRIX<T,2>& S) const
{
    MATRIX<T,3,2> U;MATRIX<T,2> V;DIAGONAL_MATRIX<T,2> D;Fast_Singular_Value_Decomposition(U,D,V);
    Q=U.Times_Transpose(V);S=SYMMETRIC_MATRIX<T,2>::Conjugate(V,D);
}
//#####################################################################
template void MATRIX<float,3,2>::Fast_Singular_Value_Decomposition(MATRIX<float,3,2>&,DIAGONAL_MATRIX<float,2>&,MATRIX<float,2>&) const;
template void MATRIX<double,3,2>::Fast_Singular_Value_Decomposition(MATRIX<double,3,2>&,DIAGONAL_MATRIX<double,2>&,MATRIX<double,2>&) const;
template void MATRIX<float,3,2>::Fast_Indefinite_Polar_Decomposition(MATRIX<float,3,2>& Q,SYMMETRIC_MATRIX<float,2>& S) const;
template void MATRIX<double,3,2>::Fast_Indefinite_Polar_Decomposition(MATRIX<double,3,2>& Q,SYMMETRIC_MATRIX<double,2>& S) const;
