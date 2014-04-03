//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/UNIFORM_ARRAY_ITERATOR.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <iomanip>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OCTAVE_OUTPUT<T>::
OCTAVE_OUTPUT(const char* file)
    :out(file)
{
    out<<std::setprecision(sizeof(T)*2);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OCTAVE_OUTPUT<T>::
~OCTAVE_OUTPUT()
{
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<class T2,class T_VECTOR> void OCTAVE_OUTPUT<T>::
Write(const char* name,const VECTOR_BASE<T2,T_VECTOR>& v)
{
    out<<"# name: "<<name<<"\n# type: matrix\n# rows: "<<v.Size()<<"\n# columns: "<<sizeof(v(1))/sizeof(typename SCALAR_POLICY<T2>::TYPE)<<"\n";
    for(int i=1;i<=v.Size();i++)
        Write_Entry(v(i));
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<class T2,class T_MATRIX> void OCTAVE_OUTPUT<T>::
Write(const char* name,const MATRIX_BASE<T2,T_MATRIX>& m)
{
    out<<"# name: "<<name<<"\n# type: matrix\n# rows: "<<m.Rows()<<"\n# columns: "<<m.Columns()<<"\n";
    for(int i=1;i<=m.Rows();i++){
        for(int j=1;j<=m.Columns();j++)
            out<<m(i,j)<<" ";
        out<<std::endl;}
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Write(const char* name,const KRYLOV_SYSTEM_BASE<T>& m,KRYLOV_VECTOR_BASE<T>& l,KRYLOV_VECTOR_BASE<T>& r)
{
    int b=r.Raw_Size();
    Begin_Sparse_Matrix(name,l.Raw_Size(),b);
    l*=0;
    r*=0;
    for(int i=1;i<=b;i++){
        r.Raw_Get(i)=1;
        m.Multiply(r,l);
        r.Raw_Get(i)=0;
        Append_Sparse_Column(l);}

    End_Sparse_Matrix();
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Write(const char* name,const KRYLOV_VECTOR_BASE<T>& v)
{
    Begin_Sparse_Matrix(name,v.Raw_Size(),1);
    Append_Sparse_Column(v);
    End_Sparse_Matrix();
}
//#####################################################################
// Function Write_Projection
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Write_Projection(const char* name,const KRYLOV_SYSTEM_BASE<T>& m,KRYLOV_VECTOR_BASE<T>& r)
{
    int b=r.Raw_Size();
    Begin_Sparse_Matrix(name,r.Raw_Size(),b);
    r*=0;
    for(int i=1;i<=b;i++){
        r*=0;
        r.Raw_Get(i)=1;
        m.Project(r);
        Append_Sparse_Column(r);}

    End_Sparse_Matrix();
}
//#####################################################################
// Function Write_Preconditioner
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Write_Preconditioner(const char* name,const KRYLOV_SYSTEM_BASE<T>& m,KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& s)
{
    int b=r.Raw_Size();
    Begin_Sparse_Matrix(name,r.Raw_Size(),b);
    r*=0;
    for(int i=1;i<=b;i++){
        r*=0;
        r.Raw_Get(i)=1;
        Append_Sparse_Column(m.Precondition(r,s));}

    End_Sparse_Matrix();
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Write(const char* name,const SPARSE_MATRIX_FLAT_MXN<T>& m)
{
    SPARSE_MATRIX_FLAT_MXN<T> t;
    m.Transpose(t);
    Write_Transpose(name,t);
}
//#####################################################################
// Function Write_Transpose
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Write_Transpose(const char* name,const SPARSE_MATRIX_FLAT_MXN<T>& m)
{
    out<<"# name: "<<name<<"\n# type: sparse matrix\n# nnz: "<<m.A.m<<"\n# rows: "<<m.n<<"\n# columns: "<<m.m<<"\n";
    for(int i=1;i<=m.m;i++){
        int s=m.offsets(i),e=m.offsets(i+1);
        for(int j=s;j<e;j++) out<<m.A(j).j<<" "<<i<<" "<<m.A(j).a<<"\n";}
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Write(const char* name,const SPARSE_MATRIX_FLAT_NXN<T>& m)
{
    SPARSE_MATRIX_FLAT_NXN<T> t;
    m.Transpose(t);
    Write_Transpose(name,t);
}
//#####################################################################
// Function Write_Transpose
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Write_Transpose(const char* name,const SPARSE_MATRIX_FLAT_NXN<T>& m)
{
    out<<"# name: "<<name<<"\n# type: sparse matrix\n# nnz: "<<m.A.m<<"\n# rows: "<<m.n<<"\n# columns: "<<m.n<<"\n";
    for(int i=1;i<=m.n;i++){
        int s=m.offsets(i),e=m.offsets(i+1);
        for(int j=s;j<e;j++) out<<(m.A(j).j)<<" "<<i<<" "<<m.A(j).a<<std::endl;}
}
//#####################################################################
// Function Begin_Sparse_Matrix
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Begin_Sparse_Matrix(const char* name,int m,int n)
{
    out<<"# name: "<<name<<"\n# type: sparse matrix\n# nnz: ";
    nnz_pos=out.tellp();
    out<<"                    \n# rows: "<<m<<"\n# columns: "<<n<<"\n";
    nnz=0;
    current_column=0;
}
//#####################################################################
// Function Add_Sparse_Entry
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Add_Sparse_Entry(int r,int c,T x)
{
    OCTAVE_SPARSE_MATRIX_ENTRY<T> e;
    e.r=r;
    e.c=c;
    e.x=x;
    internal.Append(e);
}
//#####################################################################
// Function End_Sparse_Matrix
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
End_Sparse_Matrix()
{
    if(internal.m){
        Sort(internal);
        for(int i=1;i<=internal.m;i++)
            out<<internal(i).r<<" "<<internal(i).c<<" "<<internal(i).x<<"\n";
        nnz=internal.m;
        internal.Remove_All();}

    std::streampos end=out.tellp();
    out.seekp(nnz_pos);
    out<<nnz;
    out.seekp(end);
}
//#####################################################################
// Function Append_Sparse_Column
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Append_Sparse_Column(const KRYLOV_VECTOR_BASE<T>& v)
{
    int n=v.Raw_Size();
    current_column++;
    for(int j=1;j<=n;j++)
        if(T x=v.Raw_Get(j)){
            out<<j<<" "<<current_column<<" "<<x<<std::endl;
            nnz++;}
}
//#####################################################################
// Function Append_Sparse_Column
//#####################################################################
template<class T> template<class T2,class T_ARRAY> void OCTAVE_OUTPUT<T>::
Append_Sparse_Column(const ARRAY_BASE<T2,T_ARRAY>& v)
{
    current_column++;
    for(int j=1;j<=v.Size();j++)
        if(T x=v(j)){
            out<<j<<" "<<current_column<<" "<<x<<std::endl;
            nnz++;}
}
//#####################################################################
// Function Skip_Sparse_Column
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Skip_Sparse_Column(int n)
{
    current_column+=n;
}
//#####################################################################
// Function Append_Sparse_Diagonal_Block
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Append_Sparse_Diagonal_Block(T x)
{
    current_column++;
    out<<current_column<<" "<<current_column<<" "<<x<<std::endl;
    nnz++;
}
//#####################################################################
// Function Append_Sparse_Diagonal_Block
//#####################################################################
template<class T> template<class T2,class T_MATRIX> void OCTAVE_OUTPUT<T>::
Append_Sparse_Diagonal_Block(const MATRIX_BASE<T2,T_MATRIX>& m)
{
    assert(m.Rows()==m.Columns());
    for(int i=1;i<=m.Rows();i++)
        for(int j=1;j<=m.Columns();j++)
            if(T x=m(i,j)){
                out<<(i+current_column)<<" "<<(j+current_column)<<" "<<x<<std::endl;
                nnz++;}
    current_column+=m.Rows();
}
//#####################################################################
// Function Append_Sparse_Diagonal_Block
//#####################################################################
template<class T> template<int d> void OCTAVE_OUTPUT<T>::
Append_Sparse_Diagonal_Block(const DIAGONAL_MATRIX<T,d>& m)
{
    for(int i=1;i<=d;i++){
        current_column++;
        if(T x=m(i,i)){
            out<<current_column<<" "<<current_column<<" "<<x<<std::endl;
            nnz++;}}
}
//#####################################################################
// Function Append_Sparse_Diagonal_Block
//#####################################################################
template<class T> template<int d> void OCTAVE_OUTPUT<T>::
Append_Sparse_Diagonal_Block(const SYMMETRIC_MATRIX<T,d>& m)
{
    for(int i=1;i<=d;i++)
        for(int j=1;j<=d;j++)
            if(T x=m(i,j)){
                out<<(i+current_column)<<" "<<(j+current_column)<<" "<<x<<std::endl;
                nnz++;}
    current_column+=m.Rows();
}
//#####################################################################
// Operator <
//#####################################################################
template<class T> bool OCTAVE_SPARSE_MATRIX_ENTRY<T>::
operator<(const OCTAVE_SPARSE_MATRIX_ENTRY<T>& o) const
{
    if(c!=o.c) return c<o.c;
    return r<o.r;
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<class T2> void OCTAVE_OUTPUT<T>::
Write(const char* name,const ARRAY<T2,VECTOR<int,2> >& m)
{
    out<<"# name: "<<name<<"\n# type: matrix\n# rows: "<<m.counts.y<<"\n# columns: "<<m.counts.x<<"\n";
    for(UNIFORM_ARRAY_ITERATOR<2> it(m.domain);it.Valid();it.Next()){
        out<<m(it.Index())<<" ";
        if(it.index(2)>=it.domain.max_corner(2)) out<<"\n";}
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<class T2,int d> void OCTAVE_OUTPUT<T>::
Write(const char* name,const ARRAY<VECTOR<T2,d>,VECTOR<int,2> >& m)
{
    out<<"# name: "<<name<<"\n# type: matrix\n# ndims: 3\n"<<m.counts.x<<" "<<m.counts.y<<" "<<d<<"\n";
    for(int i=1;i<=d;i++)
        for(UNIFORM_ARRAY_ITERATOR<2> it(m.domain);it.Valid();it.Next())
            out<<m(it.Index())(i)<<std::endl;
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<class T2,class T_ARRAY> void OCTAVE_OUTPUT<T>::
Write(const char* name,const ARRAY_BASE<T2,T_ARRAY>& v,int)
{
    out<<"# name: "<<name<<"\n# type: matrix\n# rows: "<<v.Size()<<"\n# columns: "<<sizeof(v(1))/sizeof(typename SCALAR_POLICY<T2>::TYPE)<<"\n";
    for(int i=1;i<=v.Size();i++)
        Write_Entry(v(i));
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Write(const char* name,T scalar)
{
    out<<"# name: "<<name<<"\n# type: scalar\n"<<scalar<<"\n";
}
//#####################################################################
// Function Write_Entry
//#####################################################################
template<class T> void OCTAVE_OUTPUT<T>::
Write_Entry(T x)
{
    out<<x<<"\n";
}
//#####################################################################
// Function Write_Entry
//#####################################################################
template<class T> template<class T2,int d> void OCTAVE_OUTPUT<T>::
Write_Entry(const VECTOR<T2,d>& x)
{
    Write_Raw(out,x);
    out<<"\n";
}
template class OCTAVE_OUTPUT<float>;
template void OCTAVE_OUTPUT<float>::Append_Sparse_Column<float,ARRAY_VIEW<float const,int> >(ARRAY_BASE<float,ARRAY_VIEW<float const,int>,int> const&);
template void OCTAVE_OUTPUT<float>::Append_Sparse_Column<float,ARRAY_VIEW<float,int> >(ARRAY_BASE<float,ARRAY_VIEW<float,int>,int> const&);
template void OCTAVE_OUTPUT<float>::Append_Sparse_Diagonal_Block<3>(SYMMETRIC_MATRIX<float,3> const&);
template void OCTAVE_OUTPUT<float>::Append_Sparse_Diagonal_Block<float,MATRIX<float,0,0> >(MATRIX_BASE<float,MATRIX<float,0,0> > const&);
template void OCTAVE_OUTPUT<float>::Append_Sparse_Diagonal_Block<float,MATRIX<float,1,1> >(MATRIX_BASE<float,MATRIX<float,1,1> > const&);
template void OCTAVE_OUTPUT<float>::Write<VECTOR<int,1>,ARRAY_VIEW<VECTOR<int,1> const,int> >(char const*,ARRAY_BASE<VECTOR<int,1>,ARRAY_VIEW<VECTOR<int,1> const,int>,int> const&,int);
template void OCTAVE_OUTPUT<float>::Write<VECTOR<int,2>,ARRAY_VIEW<VECTOR<int,2> const,int> >(char const*,ARRAY_BASE<VECTOR<int,2>,ARRAY_VIEW<VECTOR<int,2> const,int>,int> const&,int);
template void OCTAVE_OUTPUT<float>::Write<VECTOR<int,3>,ARRAY_VIEW<VECTOR<int,3> const,int> >(char const*,ARRAY_BASE<VECTOR<int,3>,ARRAY_VIEW<VECTOR<int,3> const,int>,int> const&,int);
template void OCTAVE_OUTPUT<float>::Write<VECTOR<int,4>,ARRAY_VIEW<VECTOR<int,4> const,int> >(char const*,ARRAY_BASE<VECTOR<int,4>,ARRAY_VIEW<VECTOR<int,4> const,int>,int> const&,int);
template void OCTAVE_OUTPUT<float>::Write<float,ARRAY_VIEW<float const,int> >(char const*,ARRAY_BASE<float,ARRAY_VIEW<float const,int>,int> const&,int);
template void OCTAVE_OUTPUT<float>::Write<float,VECTOR_ND<float> >(char const*,VECTOR_BASE<float,VECTOR_ND<float> > const&);
template void OCTAVE_OUTPUT<float>::Write<float>(char const*,ARRAY<float,VECTOR<int,2> > const&);
template void OCTAVE_OUTPUT<float>::Write<int,ARRAY_VIEW<int const,int> >(char const*,ARRAY_BASE<int,ARRAY_VIEW<int const,int>,int> const&,int);
template void OCTAVE_OUTPUT<float>::Write<float,ARRAY_VIEW<float,int> >(char const*,ARRAY_BASE<float,ARRAY_VIEW<float,int>,int> const&,int);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OCTAVE_OUTPUT<double>;
template void OCTAVE_OUTPUT<double>::Append_Sparse_Column<double,ARRAY_VIEW<double const,int> >(ARRAY_BASE<double,ARRAY_VIEW<double const,int>,int> const&);
template void OCTAVE_OUTPUT<double>::Append_Sparse_Column<double,ARRAY_VIEW<double,int> >(ARRAY_BASE<double,ARRAY_VIEW<double,int>,int> const&);
template void OCTAVE_OUTPUT<double>::Append_Sparse_Diagonal_Block<3>(SYMMETRIC_MATRIX<double,3> const&);
template void OCTAVE_OUTPUT<double>::Append_Sparse_Diagonal_Block<double,MATRIX<double,0,0> >(MATRIX_BASE<double,MATRIX<double,0,0> > const&);
template void OCTAVE_OUTPUT<double>::Append_Sparse_Diagonal_Block<double,MATRIX<double,1,1> >(MATRIX_BASE<double,MATRIX<double,1,1> > const&);
template void OCTAVE_OUTPUT<double>::Write<VECTOR<int,1>,ARRAY_VIEW<VECTOR<int,1> const,int> >(char const*,ARRAY_BASE<VECTOR<int,1>,ARRAY_VIEW<VECTOR<int,1> const,int>,int> const&,int);
template void OCTAVE_OUTPUT<double>::Write<VECTOR<int,2>,ARRAY_VIEW<VECTOR<int,2> const,int> >(char const*,ARRAY_BASE<VECTOR<int,2>,ARRAY_VIEW<VECTOR<int,2> const,int>,int> const&,int);
template void OCTAVE_OUTPUT<double>::Write<VECTOR<int,3>,ARRAY_VIEW<VECTOR<int,3> const,int> >(char const*,ARRAY_BASE<VECTOR<int,3>,ARRAY_VIEW<VECTOR<int,3> const,int>,int> const&,int);
template void OCTAVE_OUTPUT<double>::Write<VECTOR<int,4>,ARRAY_VIEW<VECTOR<int,4> const,int> >(char const*,ARRAY_BASE<VECTOR<int,4>,ARRAY_VIEW<VECTOR<int,4> const,int>,int> const&,int);
template void OCTAVE_OUTPUT<double>::Write<double,3>(char const*,ARRAY<VECTOR<double,3>,VECTOR<int,2> > const&);
template void OCTAVE_OUTPUT<double>::Write<double,ARRAY_VIEW<double const,int> >(char const*,ARRAY_BASE<double,ARRAY_VIEW<double const,int>,int> const&,int);
template void OCTAVE_OUTPUT<double>::Write<double,VECTOR_ND<double> >(char const*,VECTOR_BASE<double,VECTOR_ND<double> > const&);
template void OCTAVE_OUTPUT<double>::Write<int,ARRAY_VIEW<int const,int> >(char const*,ARRAY_BASE<int,ARRAY_VIEW<int const,int>,int> const&,int);
template void OCTAVE_OUTPUT<double>::Write<double,ARRAY_VIEW<double,int> >(char const*,ARRAY_BASE<double,ARRAY_VIEW<double,int>,int> const&,int);
template void OCTAVE_OUTPUT<double>::Write<double>(char const*,ARRAY<double,VECTOR<int,2> > const&);
#endif
