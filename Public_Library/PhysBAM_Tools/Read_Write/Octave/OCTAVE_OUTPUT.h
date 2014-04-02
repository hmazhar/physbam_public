//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace OCTAVE_OUTPUT 
//#####################################################################
#ifndef __OCTAVE_OUTPUT__
#define __OCTAVE_OUTPUT__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <fstream>
namespace PhysBAM{
template<class T,class T_VECTOR> class VECTOR_BASE;
template<class T,class T_MATRIX> class MATRIX_BASE;
template<class T,class T_ARRAY,class ID> class ARRAY_BASE;
template<class T,int d> class SYMMETRIC_MATRIX;
template<class T,int d> class DIAGONAL_MATRIX;
template<class T> class KRYLOV_SYSTEM_BASE;
template<class T> class KRYLOV_VECTOR_BASE;
template<class T> class SPARSE_MATRIX_FLAT_MXN;
template<class T> class SPARSE_MATRIX_FLAT_NXN;
template<class T>
struct OCTAVE_SPARSE_MATRIX_ENTRY
{
    int r,c;
    T x;
    bool operator<(const OCTAVE_SPARSE_MATRIX_ENTRY<T>& o) const;
};

template<class T>
class OCTAVE_OUTPUT{
public:
    OCTAVE_OUTPUT(const char* file);
    ~OCTAVE_OUTPUT();

    std::ofstream out;
private:
    ARRAY<OCTAVE_SPARSE_MATRIX_ENTRY<T> > internal;
    std::streampos nnz_pos;
    int nnz;
    int current_column;

public:

    template<class T2,class T_ARRAY,class ID> void Append_Sparse_Column(const ARRAY_BASE<T2,T_ARRAY,ID>& v)
    {Append_Sparse_Column(v.Array_View(ID(1),Value(v.Size())));}

    template<class T2,class ID> void Write(const char* name,ARRAY_VIEW<T2,ID> v)
    {Write(name,ARRAY_VIEW<T2>(Value(v.Size()),v.Get_Array_Pointer()),0);}

    template<class T2,class ID> void Write(const char* name,const ARRAY<T2,ID>& v)
    {Write(name,ARRAY_VIEW<const T2>(Value(v.Size()),v.Get_Array_Pointer()),0);}

    void Write_Entry(T x);
    template<class T2,int d> void Write_Entry(const VECTOR<T2,d>& x);
    template<class T2,class T_ARRAY> void Write(const char* name,const ARRAY_BASE<T2,T_ARRAY>& v,int helper=0);
    template<class T2,class T_VECTOR> void Write(const char* name,const VECTOR_BASE<T2,T_VECTOR>& v);
    template<class T2,class T_MATRIX> void Write(const char* name,const MATRIX_BASE<T2,T_MATRIX>& m);
    void Write(const char* name,T scalar);
    void Write(const char* name,const KRYLOV_SYSTEM_BASE<T>& m,KRYLOV_VECTOR_BASE<T>& l,KRYLOV_VECTOR_BASE<T>& r);
    void Write(const char* name,const KRYLOV_VECTOR_BASE<T>& v);
    void Write_Projection(const char* name,const KRYLOV_SYSTEM_BASE<T>& m,KRYLOV_VECTOR_BASE<T>& r);
    void Write_Preconditioner(const char* name,const KRYLOV_SYSTEM_BASE<T>& m,KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& s);
    void Write(const char* name,const SPARSE_MATRIX_FLAT_MXN<T>& m);
    void Write(const char* name,const SPARSE_MATRIX_FLAT_NXN<T>& m);
    void Write_Transpose(const char* name,const SPARSE_MATRIX_FLAT_MXN<T>& m);
    void Write_Transpose(const char* name,const SPARSE_MATRIX_FLAT_NXN<T>& m);
    void Begin_Sparse_Matrix(const char* name,int m,int n);
    void Add_Sparse_Entry(int r,int c,T x);
    void Append_Sparse_Column(const KRYLOV_VECTOR_BASE<T>& v);
    template<class T2,class T_ARRAY> void Append_Sparse_Column(const ARRAY_BASE<T2,T_ARRAY>& v);
    void Append_Sparse_Diagonal_Block(T x);
    template<class T2,class T_MATRIX> void Append_Sparse_Diagonal_Block(const MATRIX_BASE<T2,T_MATRIX>& m);
    template<int d> void Append_Sparse_Diagonal_Block(const DIAGONAL_MATRIX<T,d>& m);
    template<int d> void Append_Sparse_Diagonal_Block(const SYMMETRIC_MATRIX<T,d>& m);
    void Skip_Sparse_Column(int n=1);
    void End_Sparse_Matrix();
    template<class T2> void Write(const char* name,const ARRAY<T2,VECTOR<int,2> >& m);
    template<class T2,int d> void Write(const char* name,const ARRAY<VECTOR<T2,d>,VECTOR<int,2> >& m);
};
}
#endif
