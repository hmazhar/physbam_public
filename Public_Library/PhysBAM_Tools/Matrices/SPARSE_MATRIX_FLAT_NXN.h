//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Frederic Gibou, Jon Gretarsson, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_MATRIX_FLAT_NXN
//#####################################################################
#ifndef __SPARSE_MATRIX_FLAT_NXN__
#define __SPARSE_MATRIX_FLAT_NXN__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
namespace PhysBAM{

template<class T>
struct SPARSE_MATRIX_ENTRY
{
    int j;T a;
    SPARSE_MATRIX_ENTRY():j(0),a(0){}
    SPARSE_MATRIX_ENTRY(int index,T value):j(index),a(value){}
    bool operator<(const SPARSE_MATRIX_ENTRY& s) const {return j<s.j;}
};

template<class T>
class SPARSE_MATRIX_FLAT_NXN
{
public:
    typedef T SCALAR;
    int n; // size of the n by n matrix
    ARRAY<int> offsets;
    ARRAY<SPARSE_MATRIX_ENTRY<T> > A;
    ARRAY<int> diagonal_index;
    SPARSE_MATRIX_FLAT_NXN<T>* C;
    THREAD_QUEUE* thread_queue;

    void Set_Element(const int i,const int j,const T element)
    {(*this)(i,j)=element;}

    void Add_Element(const int i,const int j,const T element)
    {(*this)(i,j)+=element;}

    void Set_Symmetric_Elements(const int i,const int j,const T element)
    {assert(i!=j);Set_Element(i,j,element);Set_Element(j,i,element);}

    void Add_Symmetric_Elements(const int i,const int j,const T element)
    {assert(i!=j);Add_Element(i,j,element);Add_Element(j,i,element);}

    const T operator()(const int i,const int j) const
    {int index=Find_Index(i,j);assert(A(index).j==j);return A(index).a;}

//#####################################################################
    SPARSE_MATRIX_FLAT_NXN();
    SPARSE_MATRIX_FLAT_NXN(const SPARSE_MATRIX_FLAT_NXN<T>& matrix);
    SPARSE_MATRIX_FLAT_NXN(const ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& matrices);
    ~SPARSE_MATRIX_FLAT_NXN();
    SPARSE_MATRIX_FLAT_NXN& operator=(const SPARSE_MATRIX_FLAT_NXN& matrix);
    SPARSE_MATRIX_FLAT_NXN<T>* Create_Submatrix(const INTERVAL<int>& rows);
    void Set_Row_Lengths(const ARRAY<int>& lengths);
    int Find_Index(const int i,const int j) const;
    T& operator()(const int i,const int j);
    bool Element_Present(const int i,const int j) const;
    void Initialize_Diagonal_Index();
    void Times(const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const;
    void Times(const INTERVAL<int>& interval,const ARRAY<INTERVAL<int> >& ghost_intervals,const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const;
    void Times(const int row_start,const int row_end,const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const;
    void Times_Threaded(const VECTOR_ND<T>& x,VECTOR_ND<T>& result,int row_start,int row_end);
    void Negate();
    SPARSE_MATRIX_FLAT_NXN<T>& operator*=(const T a);
    SPARSE_MATRIX_FLAT_NXN<T>& operator+=(const T a);
    SPARSE_MATRIX_FLAT_NXN<T>& operator-=(const T a);
    bool Symmetric(const T tolerance=1e-7) const;
    bool Positive_Diagonal_And_Nonnegative_Row_Sum(const T tolerance=1e-7) const;
    void Transpose(SPARSE_MATRIX_FLAT_NXN<T>& A_transpose) const;
    bool Is_Transpose(const SPARSE_MATRIX_FLAT_NXN<T>& A_transpose,const T tolerance) const;
    void Solve_Forward_Substitution(const VECTOR_ND<T>& b,VECTOR_ND<T>& x,const bool diagonal_is_identity=false,const bool diagonal_is_inverted=false) const;
    void Solve_Backward_Substitution(const VECTOR_ND<T>& b,VECTOR_ND<T>& x,const bool diagonal_is_identity=false,const bool diagonal_is_inverted=false) const;
    // actually an LU saving square roots, with an inverted diagonal saving divides
    void Construct_Incomplete_Cholesky_Factorization(const bool modified_version=true,const T modified_coefficient=.97,const T zero_tolerance=1e-8,const T zero_replacement=1e-8);
    // actually an LU saving square roots, with an inverted diagonal saving divides
    void In_Place_Incomplete_Cholesky_Factorization(const bool modified_version=true,const T modified_coefficient=.97,const T zero_tolerance=1e-8,const T zero_replacement=1e-8);
    void Gauss_Seidel_Single_Iteration(VECTOR_ND<T>& x,const VECTOR_ND<T>& b);
    void Gauss_Seidel_Solve(VECTOR_ND<T>& x,const VECTOR_ND<T>& b,const T tolerance=1e-12,const int max_iterations=1000000);
    void Write_Row_Lengths();
    void Print_Row(const int row);
    void Reset();
    void Append_Entry_To_Current_Row(const int c,const T a);
    void Finish_Row();
    void Sort_Entries();
    void Conjugate_With_Diagonal_Matrix(VECTOR_ND<T>& x);
//#####################################################################
};
template<class T> inline std::ostream& operator<<(std::ostream& output_stream,const SPARSE_MATRIX_FLAT_NXN<T>& A)
{for(int i=1;i<=A.n;i++){
    for(int j=1;j<=A.n;j++)output_stream<<(A.Element_Present(i,j)?A(i,j):0)<<" ";
    output_stream<<std::endl;}
return output_stream;}
}
#endif
