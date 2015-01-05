//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Frederic Gibou, Jon Gretarsson, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_MATRIX_FLAT_NXN
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/INT_ITERATOR_THREADED.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_NXN<T>::
SPARSE_MATRIX_FLAT_NXN()
    :n(0),C(0),thread_queue(0)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_NXN<T>::
SPARSE_MATRIX_FLAT_NXN(const SPARSE_MATRIX_FLAT_NXN<T>& matrix)
    :n(matrix.n),offsets(matrix.offsets),A(matrix.A),C(0),thread_queue(0)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_NXN<T>::
SPARSE_MATRIX_FLAT_NXN(const ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& matrices)
    :n(0),C(0),thread_queue(0)
{
    for(int i=1;i<=matrices.m;i++) n+=matrices(i).n;
    offsets.Resize(n+1,false,false);offsets(1)=1;
    ARRAY<int> matrix_row_offsets(matrices.m);ARRAY<int> matrix_element_offsets(matrices.m);
    if(matrices.m){matrix_row_offsets(1)=0;matrix_element_offsets(1)=0;}
    for(int i=1;i<matrices.m;i++){matrix_row_offsets(i+1)=matrix_row_offsets(i)+matrices(i).n;matrix_element_offsets(i+1)=matrix_element_offsets(i)+matrices(i).A.m;}
    int index=1;
    for(int i=1;i<=matrices.m;i++)for(int j=1;j<=matrices(i).n;j++,index++) offsets(index+1)=offsets(index)+matrices(i).offsets(j+1)-matrices(i).offsets(j);
    A.Resize(offsets(n+1)-1);
    for(int i=1;i<=matrices.m;i++)for(int j=1;j<=matrices(i).A.m;j++){A(matrix_element_offsets(i)+j).j=matrices(i).A(j).j+matrix_row_offsets(i);A(matrix_element_offsets(i)+j).a=matrices(i).A(j).a;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_NXN<T>::
~SPARSE_MATRIX_FLAT_NXN()
{
    delete C;
}
//#####################################################################
// operator=
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_NXN<T>& SPARSE_MATRIX_FLAT_NXN<T>::
operator=(const SPARSE_MATRIX_FLAT_NXN& matrix)
{
    n=matrix.n;offsets=matrix.offsets;A=matrix.A;diagonal_index=matrix.diagonal_index;delete C;C=0;
    return *this;
}
//#####################################################################
// Function Create_Submatrix
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_NXN<T>* SPARSE_MATRIX_FLAT_NXN<T>::
Create_Submatrix(const INTERVAL<int>& rows)
{
    int entries=0;for(int index=offsets(rows.min_corner);index<offsets(rows.max_corner+1);index++)if(rows.Lazy_Inside(A(index).j)) entries++;
    SPARSE_MATRIX_FLAT_NXN<T>* submatrix=new SPARSE_MATRIX_FLAT_NXN<T>();
    submatrix->n=rows.Size()+1;
    submatrix->offsets.Resize(submatrix->n+1);
    submatrix->A.Resize(entries);
    int next_index=1,shift=rows.min_corner-1;
    for(int i=1;i<=submatrix->n;i++){
        submatrix->offsets(i)=next_index;
        for(int old_index=offsets(i+shift);old_index<offsets(i+shift+1);old_index++)if(rows.Lazy_Inside(A(old_index).j)){
            submatrix->A(next_index).j=A(old_index).j-shift;submatrix->A(next_index).a=A(old_index).a;next_index++;}}
    submatrix->offsets(submatrix->n+1)=next_index;
    //std::cout<<"Offsets "<<submatrix->offsets<<std::endl;
    return submatrix;
}
//#####################################################################
// Function Set_Row_Lengths
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Set_Row_Lengths(const ARRAY<int>& lengths)
{
    diagonal_index.Clean_Memory();delete C;C=0;n=lengths.m;
    offsets.Resize(n+1,false,false);offsets(1)=1;
    for(int i=1;i<=n;i++){assert(lengths(i));offsets(i+1)=offsets(i)+lengths(i);}
    A.Resize(offsets(n+1)-1);
}
//#####################################################################
// Function Find_Index
//#####################################################################
template<class T> int SPARSE_MATRIX_FLAT_NXN<T>::
Find_Index(const int i,const int j) const
{
    assert(A.m);assert(1<=i && i<=n);assert(1<=j && j<=n);
    int index=offsets(i);while(A(index).j && A(index).j<j)index++;
    assert(index<offsets(i+1));return index;
}
//#####################################################################
// Function operator()
//#####################################################################
template<class T> T& SPARSE_MATRIX_FLAT_NXN<T>::
operator()(const int i,const int j)
{
    int index=Find_Index(i,j);
    if(A(index).j!=j){ // need to add entry
        if(A(index).j){ // shift entries over to make room
            assert(!A(offsets(i+1)-1).j);
            for(int jj=offsets(i+1)-1;jj>index;jj--)A(jj)=A(jj-1);}
        A(index).j=j;A(index).a=0;}
    return A(index).a;
}
//#####################################################################
// Function Element_Present
//#####################################################################
template<class T> bool SPARSE_MATRIX_FLAT_NXN<T>::
Element_Present(const int i,const int j) const
{
    assert(1<=i && i<=n);assert(1<=j && j<=n);
    for(int index=offsets(i);index<offsets(i+1);index++)if(A(index).j==j) return true;
    return false;
}
//#####################################################################
// Function Initialize_Diagonal_Index
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Initialize_Diagonal_Index()
{
    diagonal_index.Resize(n,false,false);
    for(int i=1;i<=n;i++){diagonal_index(i)=Find_Index(i,i);assert(A(diagonal_index(i)).j==i);}
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Times(const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const
{
    //TODO: Remove const casting
    INT_ITERATOR_THREADED_ALPHA<SPARSE_MATRIX_FLAT_NXN<T> >(1,n,thread_queue).template Run<const VECTOR_ND<T>&,VECTOR_ND<T>&>(*const_cast<SPARSE_MATRIX_FLAT_NXN<T>*>(this),&SPARSE_MATRIX_FLAT_NXN<T>::Times_Threaded,x,result);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Times(const INTERVAL<int>& interval,const ARRAY<INTERVAL<int> >& ghost_intervals,const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const
{
    Times(interval.min_corner,interval.max_corner,x,result);
    for(int i=1;i<=ghost_intervals.m;i++) Times(ghost_intervals(i).min_corner,ghost_intervals(i).max_corner,x,result);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Times(const int row_start,const int row_end,const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const
{
    int index=offsets(row_start);
    for(int i=row_start;i<=row_end;i++){
        int end=offsets(i+1);T sum=0;
        for(;index<end;index++)sum+=A(index).a*x(A(index).j);
        result(i)=sum;}
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Times_Threaded(const VECTOR_ND<T>& x,VECTOR_ND<T>& result,int row_start,int row_end)
{
    Times(row_start,row_end,x,result); //TODO: Remove this function
}
//#####################################################################
// Function Negate
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Negate()
{
    for(int index=1;index<=A.m;index++) A(index).a=-A(index).a;
}
//#####################################################################
// Function operator*=
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_NXN<T>& SPARSE_MATRIX_FLAT_NXN<T>::
operator*=(const T a)
{
    for(int index=1;index<=A.m;index++) A(index).a*=a;
    return *this;
}
//#####################################################################
// Function operator+=
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_NXN<T>& SPARSE_MATRIX_FLAT_NXN<T>::
operator+=(const T a)
{
    for(int i=1;i<=n;i++) (*this)(i,i)+=a;
    return *this;
}
//#####################################################################
// Function operator-=
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_NXN<T>& SPARSE_MATRIX_FLAT_NXN<T>::
operator-=(const T a)
{
    for(int i=1;i<=n;i++) (*this)(i,i)-=a;
    return *this;
}
//#####################################################################
// Function Symmetric
//#####################################################################
template<class T> bool SPARSE_MATRIX_FLAT_NXN<T>::
Symmetric(const T tolerance) const
{
    for(int i=1;i<=n;i++)for(int index=offsets(i);index<offsets(i+1);index++)
        if(abs(A(index).a-(*this)(A(index).j,i))>tolerance) return false;
    return true;
}
//#####################################################################
// Function Positive_Diagonal_And_Nonnegative_Row_Sum
//#####################################################################
template<class T> bool SPARSE_MATRIX_FLAT_NXN<T>::
Positive_Diagonal_And_Nonnegative_Row_Sum(const T tolerance) const
{
    bool return_value=true;
    for(int i=1;i<=n;i++){
        if((*this)(i,i)<=0){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            LOG::cout<<"diagonal entry "<<i<<" contains nonpositive element: "<<(*this)(i,i)<<std::endl;
#endif
            return false;}
        T sum=0;for(int index=offsets(i);index<offsets(i+1);index++)sum+=A(index).a;
        if(sum<-tolerance){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            LOG::cout<<"sum of row "<<i<<" is negative: "<<sum<<std::endl;
#endif
            return_value=false;}}
    return return_value;
}
//#####################################################################
// Function Transpose
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Transpose(SPARSE_MATRIX_FLAT_NXN<T>& A_transpose) const
{
    ARRAY<int> row_lengths(n);for(int index=1;index<=A.m;index++)row_lengths(A(index).j)++;
    A_transpose.Set_Row_Lengths(row_lengths);
    for(int i=1;i<=n;i++) for(int index=offsets(i);index<offsets(i+1);index++) A_transpose(A(index).j,i)=A(index).a;
}
//#####################################################################
// Function Is_Transpose
//#####################################################################
template<class T> bool SPARSE_MATRIX_FLAT_NXN<T>::
Is_Transpose(const SPARSE_MATRIX_FLAT_NXN<T>& A_transpose,const T tolerance) const
{
    assert(A_transpose.n==n);
    for(int i=1;i<=n;i++)for(int index=offsets(i);index<offsets(i+1);index++)if(abs(A(index).a-A_transpose(A(index).j,i))>tolerance) return false;
    return true;
}
//#####################################################################
// Function Solve_Forward_Substitution
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Solve_Forward_Substitution(const VECTOR_ND<T>& b,VECTOR_ND<T>& x,const bool diagonal_is_identity,const bool diagonal_is_inverted) const
{
    if(diagonal_is_identity) for(int i=1;i<=n;i++){
        T sum=0;for(int index=offsets(i);index<diagonal_index(i);index++)sum+=A(index).a*x(A(index).j);
        x(i)=b(i)-sum;}
    else if(!diagonal_is_inverted) for(int i=1;i<=n;i++){
        T sum=0;for(int index=offsets(i);index<diagonal_index(i);index++)sum+=A(index).a*x(A(index).j);
        x(i)=(b(i)-sum)/A(diagonal_index(i)).a;}
    else for(int i=1;i<=n;i++){
        T sum=0;for(int index=offsets(i);index<diagonal_index(i);index++)sum+=A(index).a*x(A(index).j);
        x(i)=(b(i)-sum)*A(diagonal_index(i)).a;}
}
//#####################################################################
// Function Solve_Backward_Substitution
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Solve_Backward_Substitution(const VECTOR_ND<T>& b,VECTOR_ND<T>& x,const bool diagonal_is_identity,const bool diagonal_is_inverted) const
{
    if(diagonal_is_identity) for(int i=n;i>=1;i--){
        T sum=0;for(int index=diagonal_index(i)+1;index<offsets(i+1);index++)sum+=A(index).a*x(A(index).j);
        x(i)=b(i)-sum;}
    else if(!diagonal_is_inverted) for(int i=n;i>=1;i--){
        T sum=0;for(int index=diagonal_index(i)+1;index<offsets(i+1);index++)sum+=A(index).a*x(A(index).j);
        x(i)=(b(i)-sum)/A(diagonal_index(i)).a;}
    else for(int i=n;i>=1;i--){
        T sum=0;for(int index=diagonal_index(i)+1;index<offsets(i+1);index++)sum+=A(index).a*x(A(index).j);
        x(i)=(b(i)-sum)*A(diagonal_index(i)).a;}
}
//#####################################################################
// Function Construct_Incomplete_Cholesky_Factorization
//#####################################################################
// actually an LU saving square roots, with an inverted diagonal saving divides
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Construct_Incomplete_Cholesky_Factorization(const bool modified_version,const T modified_coefficient,const T zero_tolerance,const T zero_replacement)
{
    delete C;C=new SPARSE_MATRIX_FLAT_NXN<T>(*this);
    C->In_Place_Incomplete_Cholesky_Factorization(modified_version,modified_coefficient,zero_tolerance,zero_replacement);
}
//#####################################################################
// Function In_Place_Incomplete_Cholesky_Factorization
//#####################################################################
// actually an LU saving square roots, with an inverted diagonal saving divides
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
In_Place_Incomplete_Cholesky_Factorization(const bool modified_version,const T modified_coefficient,const T zero_tolerance,const T zero_replacement)
{
    Initialize_Diagonal_Index();
    for(int i=1;i<=n;i++){ // for each row
        int row_diagonal_index=diagonal_index(i),row_end=offsets(i+1)-1;T sum=0;
        for(int k_bar=offsets(i);k_bar<row_diagonal_index;k_bar++){ // for all the entries before the diagonal element
            int k=A(k_bar).j;int row2_diagonal_index=diagonal_index(k),row2_end=offsets(k+1)-1;
            A(k_bar).a*=A(row2_diagonal_index).a; // divide by the diagonal element (which has already been inverted)
            int j_bar=k_bar+1; // start with the next element in the row, when subtracting the dot product
            for(int i_bar=row2_diagonal_index+1;i_bar<=row2_end;i_bar++){ // run through the rest of the elements in the row2
                int i=A(i_bar).j;T dot_product_term=A(k_bar).a*A(i_bar).a;
                while(j_bar<row_end && A(j_bar).j<i) j_bar++; // gets j_bar such that j_bar>=i
                if(A(j_bar).j==i) A(j_bar).a-=dot_product_term;else if(modified_version) sum+=dot_product_term;}}
        T denominator=A(row_diagonal_index).a-modified_coefficient*sum;
        if(i==n && denominator<=zero_tolerance) A(row_diagonal_index).a=1/zero_replacement; // ensure last diagonal element is not zero
        else A(row_diagonal_index).a=1/denominator;} // finally, store the diagonal element in inverted form
}
//#####################################################################
// Function Gauss_Seidel_Single_Iteration
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Gauss_Seidel_Single_Iteration(VECTOR_ND<T>& x,const VECTOR_ND<T>& b)
{
    assert(x.n==b.n && x.n==n);
    for(int i=1;i<=n;i++){
        T rho=0;T diagonal_entry=0;
        for(int index=offsets(i);index<offsets(i+1);index++){
            if(A(index).j==i) diagonal_entry=A(index).a;
            else rho+=A(index).a*x(A(index).j);}
        x(i)=(b(i)-rho)/diagonal_entry;}
}
//#####################################################################
// Function Gauss_Seidel_Solve
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Gauss_Seidel_Solve(VECTOR_ND<T>& x,const VECTOR_ND<T>& b,const T tolerance,const int max_iterations)
{
    assert(x.n==b.n && x.n==n);
    VECTOR_ND<T> last_x(x);
    for(int k=1;k<=max_iterations;k++){
        Gauss_Seidel_Single_Iteration(x,b);
        T residual=0;for(int j=1;j<=n;j++){residual+=sqr(last_x(j)-x(j));last_x(j)=x(j);}if(residual < tolerance) return;}
}
//#####################################################################
// Function Write_Row_Lengths
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Write_Row_Lengths()
{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    for(int i=1;i<=n;i++)LOG::cout<<offsets(i+1)-offsets(i)<<" ";LOG::cout<<std::endl;
#endif
}
//#####################################################################
// Function Print_Row
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Print_Row(const int row)
{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    for(int i=1;i<=n;i++)if(Element_Present(row,i))LOG::cout<<"Col: "<<i<<" Val: "<<(*this)(row,i)<<",  ";
    LOG::cout<<std::endl;
#endif
}
//#####################################################################
// Function Reset
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Reset()
{
    delete C;
    C=0;
    offsets.Remove_All();
    A.Remove_All();
    offsets.Append(1);
    diagonal_index.Remove_All();
}
//#####################################################################
// Function Append_Entry_To_Current_Row
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Append_Entry_To_Current_Row(const int c,const T a)
{
    A.Append(SPARSE_MATRIX_ENTRY<T>(c,a));
}
//#####################################################################
// Function Finish_Row
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Finish_Row()
{
    offsets.Append(A.m+1);
    n++;
}
//#####################################################################
// Function Sort_Entries
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Sort_Entries()
{
    for(int i=1;i<=n;i++){ARRAY_VIEW<SPARSE_MATRIX_ENTRY<T> > view(A.Array_View(offsets(i),offsets(i+1)-offsets(i)));Sort(view);}
}
//#####################################################################
// Function Conjugate_With_Diagonal_Matrix
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_NXN<T>::
Conjugate_With_Diagonal_Matrix(VECTOR_ND<T>& x)
{
    int index=offsets(1);
    for(int i=1;i<=n;i++){
        int end=offsets(i+1);
        for(;index<end;index++) A(index).a*=x(i)*x(A(index).j);}
}
//#####################################################################
// Function operator<<
//#####################################################################
template<class T> std::ostream&
operator<<(std::ostream& output_stream,const SPARSE_MATRIX_FLAT_NXN<T>& A)
{for(int i=1;i<=A.n;i++){
    for(int j=1;j<=A.n;j++)output_stream<<(A.Element_Present(i,j)?A(i,j):0)<<" ";
    output_stream<<std::endl;}
return output_stream;}
//#####################################################################
template class SPARSE_MATRIX_FLAT_NXN<float>;
template std::ostream& operator<<(std::ostream&,const SPARSE_MATRIX_FLAT_NXN<float>&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SPARSE_MATRIX_FLAT_NXN<double>;
template std::ostream& operator<<(std::ostream&,const SPARSE_MATRIX_FLAT_NXN<double>&);
#endif
