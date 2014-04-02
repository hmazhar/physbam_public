//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Frederic Gibou, Jon Gretarsson, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_MATRIX_FLAT_MXN
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>::
SPARSE_MATRIX_FLAT_MXN()
    :m(0),n(0),Q(0),L(0)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>::
SPARSE_MATRIX_FLAT_MXN(const SPARSE_MATRIX_FLAT_MXN<T>& matrix)
    :m(matrix.m),n(matrix.n),offsets(matrix.offsets),A(matrix.A),Q(0),L(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>::
~SPARSE_MATRIX_FLAT_MXN()
{
    delete Q;delete L;
}
//#####################################################################
// Function Create_Submatrix
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_NXN<T>* SPARSE_MATRIX_FLAT_MXN<T>::
Create_Submatrix(const INTERVAL<int>& rows)
{
    assert(rows.Size()+1==m);
    int entries=0;for(int index=offsets(1);index<offsets(m+1);index++)if(rows.Lazy_Inside(A(index).j)) entries++;
    SPARSE_MATRIX_FLAT_NXN<T>* submatrix=new SPARSE_MATRIX_FLAT_NXN<T>();
    submatrix->n=rows.Size()+1;
    submatrix->offsets.Resize(submatrix->n+1);
    submatrix->A.Resize(entries);
    int next_index=1,shift=rows.min_corner-1;
    for(int i=1;i<=submatrix->n;i++){
        submatrix->offsets(i)=next_index;
        for(int old_index=offsets(i);old_index<offsets(i+1);old_index++)if(rows.Lazy_Inside(A(old_index).j)){
            submatrix->A(next_index).j=A(old_index).j-shift;submatrix->A(next_index).a=A(old_index).a;next_index++;}}
    submatrix->offsets(submatrix->n+1)=next_index;
    return submatrix;
}
//#####################################################################
// Function Set_Row_Lengths
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Set_Row_Lengths(ARRAY_VIEW<int> lengths)
{
    m=lengths.m;offsets.Resize(m+1,false,false);offsets(1)=1;
    for(int i=1;i<=m;i++){offsets(i+1)=offsets(i)+lengths(i);}
    A.Resize(offsets(m+1)-1);
}
//#####################################################################
// Function Find_Index
//#####################################################################
template<class T> int SPARSE_MATRIX_FLAT_MXN<T>::
Find_Index(const int i,const int j) const
{
    assert(A.m);assert(1<=i && i<=m);assert(1<=j && j<=n);
    int index=offsets(i);while(A(index).j && A(index).j<j)index++;
    assert(index<offsets(i+1));return index;
}
//#####################################################################
// Function Find_Index_Return
//#####################################################################
template<class T> int SPARSE_MATRIX_FLAT_MXN<T>::
Find_Index_Exists(const int i,const int j) const
{
    assert(A.m);assert(1<=i && i<=m);assert(1<=j && j<=n);
    int index=offsets(i);while(A(index).j && A(index).j<j)index++;
    if(index<offsets(i+1) && A(index).j==j)
        return index;
    else
        return 0;
}
//#####################################################################
// Function operator()
//#####################################################################
template<class T> T& SPARSE_MATRIX_FLAT_MXN<T>::
operator()(const int i,const int j)
{
    int index=Find_Index(i,j);
    if(A(index).j!=j){ // need to add entry
        if(A(index).j){ // shift entries over to make room
            assert(!A(offsets(i+1)-1).j);
            for(int jj=offsets(i+1)-1;jj>index;jj--) A(jj)=A(jj-1);}
        A(index).j=j;A(index).a=0;}
    return A(index).a;
}
//#####################################################################
// Function Element_Present
//#####################################################################
template<class T> bool SPARSE_MATRIX_FLAT_MXN<T>::
Element_Present(const int i,const int j) const
{
    assert(1<=i && i<=m);assert(1<=j && j<=n);
    for(int index=offsets(i);index<offsets(i+1);index++)if(A(index).j==j) return true;
    return false;
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Times_Add(const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const
{
    int index=offsets(1);
    for(int i=1;i<=m;i++){
        int end=offsets(i+1);T sum=(T)0;
        for(;index<end;index++) sum+=A(index).a*x(A(index).j);
        result(i)+=sum;}
}
//#####################################################################
// Function Times_Subtract
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Times_Subtract(const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const
{
    int index=offsets(1);
    for(int i=1;i<=m;i++){
        int end=offsets(i+1);T sum=(T)0;
        for(;index<end;index++) sum+=A(index).a*x(A(index).j);
        result(i)-=sum;}
}
//#####################################################################
// Function Times
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Times(const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const
{
    result.Fill(0);
    Times_Add(x,result);
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Transpose_Times_Add(const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const
{
    int index=offsets(1);
    for(int i=1;i<=m;i++){
        int end=offsets(i+1);T y=x(i);
        for(;index<end;index++) result(A(index).j)+=A(index).a*y;}
}
//#####################################################################
// Function Transpose_Times_Subtract
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Transpose_Times_Subtract(const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const
{
    int index=offsets(1);
    for(int i=1;i<=m;i++){
        int end=offsets(i+1);T y=x(i);
        for(;index<end;index++) result(A(index).j)-=A(index).a*y;}
}
//#####################################################################
// Function Times
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Transpose_Times(const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const
{
    result.Fill(0);
    Transpose_Times_Add(x,result);
}
//#####################################################################
// Function Negate
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Negate()
{
    for(int index=1;index<=A.m;index++) A(index).a=-A(index).a;
}
//#####################################################################
// Function operator*=
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>& SPARSE_MATRIX_FLAT_MXN<T>::
operator*=(const T a)
{
    for(int index=1;index<=A.m;index++) A(index).a*=a;
    return *this;
}
//#####################################################################
// Function operator/=
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>& SPARSE_MATRIX_FLAT_MXN<T>::
operator/=(const T a)
{
    return *this*=1/a;
}
//#####################################################################
// Function operator+=
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>& SPARSE_MATRIX_FLAT_MXN<T>::
operator+=(const T a)
{
    int index=offsets(1);
    for(int i=1;i<=m;i++){
        int end=offsets(i+1),found=0;
        for(;index<end;index++){if(A(index).j==i){found=1;A(index).a+=a;}}
        PHYSBAM_ASSERT(found);}
    return *this;
}
//#####################################################################
// Function operator-=
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>& SPARSE_MATRIX_FLAT_MXN<T>::
operator-=(const T a)
{
    return *this+=-a;
}
//#####################################################################
// Function Compress
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Compress(SPARSE_MATRIX_FLAT_MXN<T>& compressed)
{
    ARRAY<int> row_lengths(m);
    int index=offsets(1);
    for(int i=1;i<=m;i++){
        int end=offsets(i+1);
        while(index<end){if(!A(index).j) break;
            index++;row_lengths(i)++;}
        index=end;}
    compressed.Set_Row_Lengths(row_lengths);
    index=offsets(1);
    int compressed_index=1;
    for(int i=1;i<=m;i++){
        int end=offsets(i+1);
        while(index<end){if(!A(index).j) break;
            compressed.A(compressed_index++)=A(index);index++;}
        index=end;}
}
//#####################################################################
// Function Transpose
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Transpose(SPARSE_MATRIX_FLAT_MXN<T>& A_transpose) const
{
    ARRAY<int> row_lengths(n);for(int index=1;index<=A.m;index++) row_lengths(A(index).j)++;
    A_transpose.Set_Row_Lengths(row_lengths);A_transpose.n=m;
    for(int i=1;i<=m;i++) for(int index=offsets(i);index<offsets(i+1);index++) A_transpose(A(index).j,i)=A(index).a;
}
//#####################################################################
// Function Times_Transpose
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
Times_Transpose(const SPARSE_MATRIX_FLAT_MXN<T>& rhs)
{
    assert(rhs.n==n);
    SPARSE_MATRIX_FLAT_MXN<T> result;
    const int columns=rhs.m;const int rows=m;
    ARRAY<int> row_lengths(rows);
    ARRAY<int> right_columns_present,left_rows_present;
    for(int i=1;i<=m;i++) if(offsets(i+1)-offsets(i)) left_rows_present.Append(i);
    for(int i=1;i<=rhs.m;i++) if(rhs.offsets(i+1)-rhs.offsets(i)) right_columns_present.Append(i);
    
    for(int i=1;i<=left_rows_present.m;i++){
        int row=left_rows_present(i);
        int start=offsets(row),end=offsets(row+1);
        for(int j=1;j<=right_columns_present.m;j++){
            int col=right_columns_present(j);
            int right_index=rhs.offsets(col);int right_end=rhs.offsets(col+1);int index=start;
            while(index<end && right_index<right_end){
                if(A(index).j==rhs.A(right_index).j){row_lengths(row)++;break;}
                else if(A(index).j>rhs.A(right_index).j) right_index++;
                else index++;}}}
    result.Set_Row_Lengths(row_lengths);result.n=columns;
    for(int i=1;i<=left_rows_present.m;i++){
        int row=left_rows_present(i);
        int start=offsets(row),end=offsets(row+1);
        for(int j=1;j<=right_columns_present.m;j++){
            int col=right_columns_present(j);
            int right_index=rhs.offsets(col);int right_end=rhs.offsets(col+1);int index=start;
            while(index<end && right_index<right_end){
                if(A(index).j==rhs.A(right_index).j){result.Add_Element(row,col,A(index).a*rhs.A(right_index).a);index++;right_index++;}
                else if(A(index).j>rhs.A(right_index).j) right_index++;
                else index++;}}}
    return result;
}
//#####################################################################
// Function Times_Diagonal_Times
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
Times_Diagonal_Times(const VECTOR_ND<T> diagonal,const SPARSE_MATRIX_FLAT_MXN<T>& rhs)
{
    assert(rhs.m==n && diagonal.n==n);
    SPARSE_MATRIX_FLAT_MXN result;
    const int columns=rhs.n,rows=m;
    ARRAY<int> row_lengths(rows);
    ARRAY<int> column_indices;
    for(int row=1;row<=rows;row++){
        int start=offsets(row),end=offsets(row+1);
        column_indices.Remove_All();
        for(int i=start;i<end;i++){ // insertion sort to get actual row counts, worst-case cost C^2
            int col=A(i).j;
            for(int j=rhs.offsets(col);j<rhs.offsets(col+1);j++)
                column_indices.Append_Unique(rhs.A(j).j);
        }
        row_lengths(row)=column_indices.m;
    }

    result.Set_Row_Lengths(row_lengths);
    result.n=columns;

    for(int row=1;row<=rows;row++){
        int start=offsets(row),end=offsets(row+1);
        for(int i=start;i<end;i++){
            int col=A(i).j;
            for(int j=rhs.offsets(col);j<rhs.offsets(col+1);j++)
                result.Add_Element(row,rhs.A(j).j,A(i).a*rhs.A(j).a*diagonal(col));
        }
    }
    return result;
}
//#####################################################################
// Function Scale_Rows
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
Scale_Rows(const VECTOR_ND<T>& d) const
{
    SPARSE_MATRIX_FLAT_MXN<T> result=*this;
    for(int i=1;i<=result.m;i++)
        for(int j=result.offsets(i);j<result.offsets(i+1);j++)
            result.A(j).a*=d(i);
    return result;
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_NXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
operator+(const SPARSE_MATRIX_FLAT_NXN<T>& A_rhs) const
{
    assert(m==n && n==A_rhs.n);
    ARRAY<int> row_lengths(m);
    int left_index=offsets(1),right_index=A_rhs.offsets(1);
    for(int i=1;i<=m;i++){
        int left_end=offsets(i+1),right_end=A_rhs.offsets(i+1);
        while(left_index<left_end && right_index<right_end){
            if(A(left_index).j==A_rhs.A(right_index).j){row_lengths(i)++;left_index++;right_index++;}
            else if(A(left_index).j>A_rhs.A(right_index).j){right_index++;row_lengths(i)++;}
            else{left_index++;row_lengths(i)++;}}
        row_lengths(i)+=left_end-left_index+right_end-right_index;
        left_index=left_end;right_index=right_end;}
    SPARSE_MATRIX_FLAT_NXN<T> result;
    result.Set_Row_Lengths(row_lengths);

    int index=offsets(1);
    for(int i=1;i<=m;i++){int end=offsets(i+1);
        for(;index<end;index++) result(i,A(index).j)=A(index).a;}

    index=A_rhs.offsets(1);
    for(int i=1;i<=m;i++){int end=A_rhs.offsets(i+1);
        for(;index<end;index++) result(i,A_rhs.A(index).j)+=A_rhs.A(index).a;}
    return result;
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
operator+(const SPARSE_MATRIX_FLAT_MXN& A_rhs) const
{
    assert(m==A_rhs.m && n==A_rhs.n);
    ARRAY<int> row_lengths(m);
    int left_index=offsets(1),right_index=A_rhs.offsets(1);
    for(int i=1;i<=m;i++){
        int left_end=offsets(i+1),right_end=A_rhs.offsets(i+1);
        while(left_index<left_end && right_index<right_end){
            if(A(left_index).j==A_rhs.A(right_index).j){row_lengths(i)++;left_index++;right_index++;}
            else if(A(left_index).j>A_rhs.A(right_index).j){right_index++;row_lengths(i)++;}
            else{left_index++;row_lengths(i)++;}}
        row_lengths(i)+=left_end-left_index+right_end-right_index;
        left_index=left_end;right_index=right_end;}
    SPARSE_MATRIX_FLAT_MXN<T> result;
    result.Set_Row_Lengths(row_lengths);result.n=n;

    int index=offsets(1);
    for(int i=1;i<=m;i++){int end=offsets(i+1);
        for(;index<end;index++) result(i,A(index).j)=A(index).a;}

    index=A_rhs.offsets(1);
    for(int i=1;i<=m;i++){int end=A_rhs.offsets(i+1);
        for(;index<end;index++) result(i,A_rhs.A(index).j)+=A_rhs.A(index).a;}
    return result;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
operator-(const SPARSE_MATRIX_FLAT_MXN& A_rhs) const
{
    assert(m==A_rhs.m && n==A_rhs.n);
    ARRAY<int> row_lengths(m);
    int left_index=offsets(1),right_index=A_rhs.offsets(1);
    for(int i=1;i<=m;i++){
        int left_end=offsets(i+1),right_end=A_rhs.offsets(i+1);
        while(left_index<left_end && right_index<right_end){
            if(A(left_index).j==A_rhs.A(right_index).j){row_lengths(i)++;left_index++;right_index++;}
            else if(A(left_index).j>A_rhs.A(right_index).j){right_index++;row_lengths(i)++;}
            else{left_index++;row_lengths(i)++;}}
        row_lengths(i)+=left_end-left_index+right_end-right_index;
        left_index=left_end;right_index=right_end;}
    SPARSE_MATRIX_FLAT_MXN<T> result;
    result.Set_Row_Lengths(row_lengths);result.n=n;

    int index=offsets(1);
    for(int i=1;i<=m;i++){int end=offsets(i+1);
        for(;index<end;index++) result(i,A(index).j)=A(index).a;}

    index=A_rhs.offsets(1);
    for(int i=1;i<=m;i++){int end=A_rhs.offsets(i+1);
        for(;index<end;index++) result(i,A_rhs.A(index).j)-=A_rhs.A(index).a;}
    return result;
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
operator*(const SPARSE_MATRIX_FLAT_MXN& rhs) const
{
    assert(rhs.m==n);
    SPARSE_MATRIX_FLAT_MXN result;
    const int columns=rhs.n,rows=m;
    ARRAY<int> row_lengths(rows);
    ARRAY<int> column_indices;
    for(int row=1;row<=rows;row++){
        int start=offsets(row),end=offsets(row+1);
        column_indices.Remove_All();
        for(int i=start;i<end;i++){ // insertion sort to get actual row counts, worst-case cost C^2
            int col=A(i).j;
            for(int j=rhs.offsets(col);j<rhs.offsets(col+1);j++)
                column_indices.Append_Unique(rhs.A(j).j);
        }
        row_lengths(row)=column_indices.m;
    }

    result.Set_Row_Lengths(row_lengths);
    result.n=columns;

    for(int row=1;row<=rows;row++){
        int start=offsets(row),end=offsets(row+1);
        for(int i=start;i<end;i++){
            int col=A(i).j;
            for(int j=rhs.offsets(col);j<rhs.offsets(col+1);j++)
                result.Add_Element(row,rhs.A(j).j,A(i).a*rhs.A(j).a);
        }
    }
    return result;
}
//#####################################################################
// Function Write_Row_Lengths
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_NXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
Create_NXN_Matrix()
{
    SPARSE_MATRIX_FLAT_NXN<T> result;
    result.offsets=offsets;
    result.A=A;
    result.n=m;
    return result;
}
//#####################################################################
// Function Write_Row_Lengths
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Write_Row_Lengths()
{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    for(int i=1;i<=m;i++) LOG::cout<<offsets(i+1)-offsets(i)<<" ";
    LOG::cout<<std::endl;
#endif
}
//#####################################################################
// Function Print_Row
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Print_Row(const int row)
{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    for(int i=1;i<=n;i++) if(Element_Present(row,i)) LOG::cout<<"Col: "<<i<<" Val: "<<(*this)(row,i)<<",  ";
    LOG::cout<<std::endl;
#endif
}
//#####################################################################
// Function operator<<
//#####################################################################
template<class T> std::ostream&
operator<<(std::ostream& output_stream,const SPARSE_MATRIX_FLAT_MXN<T>& A)
{
    for(int i=1;i<=A.m;i++){
        for(int j=1;j<=A.n;j++)output_stream<<(A.Element_Present(i,j)?A(i,j):0)<<" ";
        output_stream<<std::endl;}
    return output_stream;
}
//#####################################################################
// Function Reset
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Reset(const int c)
{
    n=c;
    m=0;
    offsets.Remove_All();
    A.Remove_All();
    offsets.Append(1);
}
//#####################################################################
// Function Append_Entry_To_Current_Row
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Append_Entry_To_Current_Row(const int c,const T a)
{
    A.Append(SPARSE_MATRIX_ENTRY<T>(c,a));
}
//#####################################################################
// Function Finish_Row
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Finish_Row()
{
    offsets.Append(A.m+1);
    m++;
}
//#####################################################################
// Function Sort_Entries
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Sort_Entries()
{
    for(int i=1;i<=m;i++){ARRAY_VIEW<SPARSE_MATRIX_ENTRY<T> > view(A.Array_View(offsets(i),offsets(i+1)-offsets(i)));Sort(view);}
}
//#####################################################################
// Function Get_Row
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Get_Row(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q_i,int row)
{    
    q_i.Resize(offsets(row+1)-offsets(row));
    q_i=A.Array_View(offsets(row),offsets(row+1)-offsets(row));
}
//#####################################################################
// Function Construct_Incomplete_QR_Factorization
//#####################################################################
template<class T> void
Check_Sorted(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q)
{
    for(int i=1;i<=q.m-1;i++)
        assert(q(i).j<q(i+1).j);
}
//#####################################################################
// Function Keep_Largest_N
//#####################################################################
template<class T> void
Keep_Largest_N(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q,int n)
{
    static ARRAY<SPARSE_MATRIX_ENTRY<T> > elements;
    elements.Resize(0);
    int j;
    for(int i=1;i<=q.m;i++){
        T element_magnitude=abs(q(i).a);
        for(j=1;j<=elements.m;j++)
            if(element_magnitude>abs(elements(j).a)){
                elements.Insert(q(i),j);
                break;}
        if(j>elements.m && elements.m<n) elements.Append(q(i));
        if(elements.m>n) elements.Resize(n);}
    q=elements;
    Sort(q);
    Check_Sorted(q);
}
//#####################################################################
// Function Subtract_C_Times
//#####################################################################
template<class T> void
Subtract_C_Times(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q,T c,const ARRAY_VIEW<SPARSE_MATRIX_ENTRY<T> >& view)
{
    int q_index=0;
    int current_q_column=0;
    for(int view_index=1;view_index<=view.m;view_index++){
        while(q_index<q.m && current_q_column<view(view_index).j){q_index++;current_q_column=q(q_index).j;}
        if(current_q_column==view(view_index).j){
            q(q_index).a-=c*view(view_index).a;
            view_index++;
        }
        else if(current_q_column>view(view_index).j)
            q.Insert(SPARSE_MATRIX_ENTRY<T>(view(view_index).j,-c*view(view_index).a),q_index++);
        else
            q.Insert(SPARSE_MATRIX_ENTRY<T>(view(view_index).j,-c*view(view_index).a),q_index+1);
    }
    Check_Sorted(q);
}
//#####################################################################
// Function Two_Norm
//#####################################################################
template<class T> T
Two_Norm(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q)
{
    double total=0;
    for(int i=1;i<=q.m;i++) total+=sqr(q(i).a);
    return sqrt(total);
}
//#####################################################################
// Function Construct_Incomplete_LQ_Factorization
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Construct_Incomplete_LQ_Factorization(const int p_l,const int p_q,const T zero_tolerance,const T zero_replacement)
{
    ARRAY<SPARSE_MATRIX_ENTRY<T> > q_i,l_i;
    delete Q;delete L;Q=new SPARSE_MATRIX_FLAT_MXN<T>();L=new SPARSE_MATRIX_FLAT_MXN<T>();
    (*Q).Reset(n);(*L).Reset(m);
    for(int i=1;i<=m;i++){
        Get_Row(q_i,i);
        (*Q).Fast_Sparse_Multiply(q_i,l_i);
        Keep_Largest_N(l_i,p_l);
        for(int j=1;j<=l_i.m;j++)
            Subtract_C_Times(q_i,l_i(j).a,(*Q).A.Array_View((*Q).offsets(l_i(j).j),(*Q).offsets(l_i(j).j+1)-(*Q).offsets(l_i(j).j)));
        Keep_Largest_N(q_i,p_q);
        T magnitude=Two_Norm(q_i);
        if(magnitude<zero_tolerance){PHYSBAM_FATAL_ERROR("Zero pivot in LQ factorization");magnitude=zero_replacement;}
        T one_over_magnitude=Inverse(magnitude);
        for(int j=1;j<=q_i.m;j++)
            (*Q).Append_Entry_To_Current_Row(q_i(j).j,q_i(j).a*one_over_magnitude);
        (*Q).Finish_Row();
        for(int j=1;j<=l_i.m;j++)
            (*L).Append_Entry_To_Current_Row(l_i(j).j,l_i(j).a);
        (*L).Append_Entry_To_Current_Row(i,magnitude);
        (*L).Finish_Row();
    }
}
//#####################################################################
// Function Fast_Sparse_Multiply
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Fast_Sparse_Multiply(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q,ARRAY<SPARSE_MATRIX_ENTRY<T> >& l)
{
    SPARSE_MATRIX_FLAT_MXN<T> transpose;
    Transpose(transpose);
    l.Resize(0);
    for(int i=1;i<=q.m;i++)
        Subtract_C_Times(l,-q(i).a,transpose.A.Array_View(transpose.offsets(q(i).j),transpose.offsets(q(i).j+1)-transpose.offsets(q(i).j)));
}
//#####################################################################
template class SPARSE_MATRIX_FLAT_MXN<float>;
template std::ostream& operator<<(std::ostream&,const SPARSE_MATRIX_FLAT_MXN<float>&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SPARSE_MATRIX_FLAT_MXN<double>;
template std::ostream& operator<<(std::ostream&,const SPARSE_MATRIX_FLAT_MXN<double>&);
#endif
}
