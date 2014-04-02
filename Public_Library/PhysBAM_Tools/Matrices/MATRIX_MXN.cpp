//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Huamin Wang, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_MXN
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <limits>
using namespace PhysBAM;
//#####################################################################
// Function Jacobi_Singular_Value_Decomposition
//#####################################################################
namespace{
template<class T> static void 
Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(const MATRIX_MXN<T>& A,ARRAY<PAIR<int,T> >& max_off_diagonal_element_of_row,const int i)
{max_off_diagonal_element_of_row(i)=PAIR<int,T>();
for(int j=1;j<=A.n;j++) if(i!=j) 
    if(abs(A(i,j))>max_off_diagonal_element_of_row(i).y){max_off_diagonal_element_of_row(i).y=abs(A(i,j));max_off_diagonal_element_of_row(i).x=j;}}
template<class T> static void 
Update_Max_Off_Diagonal_Element_Of_Row_After_Column_Change(const MATRIX_MXN<T>& A,ARRAY<PAIR<int,T> >& max_off_diagonal_element_of_row,const int j)
{for(int i=1;i<=A.m;i++) if(i!=j){
    if(abs(A(i,j))>max_off_diagonal_element_of_row(i).y){max_off_diagonal_element_of_row(i).y=abs(A(i,j));max_off_diagonal_element_of_row(i).x=j;}
    else if(max_off_diagonal_element_of_row(i).x==j) Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(A,max_off_diagonal_element_of_row,i);}}
}
template<class T> void MATRIX_MXN<T>::
Jacobi_Singular_Value_Decomposition(ARRAY<VECTOR<int,2> >& left_givens_pairs,ARRAY<VECTOR<T,2> >& left_givens_coefficients,
    ARRAY<VECTOR<int,2> >& right_givens_pairs,ARRAY<VECTOR<T,2> >& right_givens_coefficients,const T tolerance,const int max_iterations)
{
    assert(m>=2 && n>=2);
    ARRAY<PAIR<int,T> > max_off_diagonal_element_of_row(m);
    for(int i=1;i<=m;i++) Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(*this,max_off_diagonal_element_of_row,i);
    left_givens_pairs.Remove_All();left_givens_coefficients.Remove_All();
    right_givens_pairs.Remove_All();right_givens_coefficients.Remove_All();
    for(int iteration=1;iteration<=max_iterations;iteration++){
        T max_off_diagonal_element=0;int i_max=0,j_max=0;
        for(int i=1;i<=m;i++)
            if(max_off_diagonal_element_of_row(i).y>max_off_diagonal_element){max_off_diagonal_element=max_off_diagonal_element_of_row(i).y;i_max=i;j_max=max_off_diagonal_element_of_row(i).x;}
        if(max_off_diagonal_element<tolerance) return;
        if(i_max>n){
            int i=j_max,j=i_max;T c,s;VECTOR<T,2>((*this)(i,i),(*this)(j,i)).Normalized().Get(c,s);
            Left_Givens_Rotation(i,j,c,-s);left_givens_pairs.Append(VECTOR<int,2>(i,j));left_givens_coefficients.Append(VECTOR<T,2>(c,s));
            Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(*this,max_off_diagonal_element_of_row,i);
            Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(*this,max_off_diagonal_element_of_row,j);}
        else if(j_max>m){
            int i=i_max,j=j_max;T c,s;VECTOR<T,2>((*this)(i,i),(*this)(i,j)).Normalized().Get(c,s);
            Right_Givens_Rotation(i,j,c,-s);right_givens_pairs.Append(VECTOR<int,2>(i,j));right_givens_coefficients.Append(VECTOR<T,2>(c,s));
            Update_Max_Off_Diagonal_Element_Of_Row_After_Column_Change(*this,max_off_diagonal_element_of_row,i);
            Update_Max_Off_Diagonal_Element_Of_Row_After_Column_Change(*this,max_off_diagonal_element_of_row,j);}
        else{
            int i=(i_max<j_max)?i_max:j_max,j=(i_max<j_max)?j_max:i_max;
            MATRIX<T,2> B((*this)(i,i),(*this)(j,i),(*this)(i,j),(*this)(j,j)),U,V;DIAGONAL_MATRIX<T,2> sigma;B.Fast_Singular_Value_Decomposition(U,sigma,V);
            T c_left,s_left,c_right,s_right;U.Column(1).Get(c_left,s_left);V.Column(1).Get(c_right,s_right);
            Left_Givens_Rotation(i,j,c_left,-s_left);left_givens_pairs.Append(VECTOR<int,2>(i,j));left_givens_coefficients.Append(VECTOR<T,2>(c_left,s_left));
            Right_Givens_Rotation(i,j,c_right,-s_right);right_givens_pairs.Append(VECTOR<int,2>(i,j));right_givens_coefficients.Append(VECTOR<T,2>(c_right,s_right));
            Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(*this,max_off_diagonal_element_of_row,i);
            Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(*this,max_off_diagonal_element_of_row,j);
            Update_Max_Off_Diagonal_Element_Of_Row_After_Column_Change(*this,max_off_diagonal_element_of_row,i);
            Update_Max_Off_Diagonal_Element_Of_Row_After_Column_Change(*this,max_off_diagonal_element_of_row,j);}}
}
//#####################################################################
template class MATRIX_MXN<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MATRIX_MXN<double>;
#endif
