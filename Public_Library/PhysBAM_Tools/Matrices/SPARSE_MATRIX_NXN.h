//#####################################################################
// Copyright 2003-2005, Ron Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_MATRIX_NXN
//#####################################################################
#ifndef __SPARSE_MATRIX_NXN__
#define __SPARSE_MATRIX_NXN__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <cassert>
namespace PhysBAM{

template<class T>
class SPARSE_MATRIX_NXN
{
public:
    int n; // size of the n by n matrix
    ARRAY<SPARSE_VECTOR_ND<T>*> A;
    VECTOR_ND<int>* diagonal_index;
    SPARSE_MATRIX_NXN<T>* C;

    SPARSE_MATRIX_NXN()
        :n(0),diagonal_index(0),C(0)
    {
        Set_Size(0);
    }

    SPARSE_MATRIX_NXN(const int n_input)
        :n(0),diagonal_index(0),C(0)
    {
        Set_Size(n_input);
    }

    SPARSE_MATRIX_NXN(const SPARSE_MATRIX_NXN& matrix)
        :n(matrix.n),diagonal_index(0),C(0)
    {
        A.Resize(n,false,false);
        for(int i=1;i<=n;i++) A(i)=new SPARSE_VECTOR_ND<T>(*matrix.A(i));
    }

    ~SPARSE_MATRIX_NXN()
    {for(int i=1;i<=n;i++) delete A(i);delete diagonal_index;delete C;}

    void Clean_Memory()
    {A.Delete_Pointers_And_Clean_Memory();delete diagonal_index;diagonal_index=0;delete C;C=0;}

    SPARSE_MATRIX_NXN& operator=(const SPARSE_MATRIX_NXN& matrix)
    {for(int i=1;i<=n;i++) delete A(i);delete diagonal_index;delete C;
    n=matrix.n;diagonal_index=0;C=0;A.Resize(n,false,false);
    for(int i=1;i<=n;i++) A(i)=new SPARSE_VECTOR_ND<T>(*matrix.A(i));
    return *this;}

    void Set_Size(const int n_input)
    {for(int i=1;i<=n;i++) delete A(i);delete diagonal_index;delete C;
    n=n_input;diagonal_index=0;C=0;A.Resize(n,false,false);
    for(int i=1;i<=n;i++) A(i)=new SPARSE_VECTOR_ND<T>(n);}

    const T operator()(const int i,const int j) const
    {assert(i>=1 && i<=n);assert(j>=1 && j<=n);return (*A(i))(j);}

    void Set_Element(const int i,const int j,const T element)
    {assert(i>=1 && i<=n);assert(j>=1 && j<=n);
    A(i)->Set_Element(j,element);}

    void Set_Symmetric_Elements(const int i,const int j,const T element)
    {assert(i!=j);Set_Element(i,j,element);Set_Element(j,i,element);}

    void Add_Element(const int i,const int j,const T element)
    {assert(i>=1 && i<=n);assert(j>=1 && j<=n);
    A(i)->Add_Element(j,element);}

    void Add_Symmetric_Elements(const int i,const int j,const T element)
    {assert(i!=j);Add_Element(i,j,element);Add_Element(j,i,element);}

    bool Element_Present(const int i,const int j)
    {assert(i>=1 && i<=n);assert(j>=1 && j<=n);
    return A(i)->Element_Present(j);}

    void Clear_Row(const int i)
    {assert(i>=1 && i<=n);A(i)->Clear();}

    void Initialize_Diagonal_Index()
    {if(!diagonal_index) diagonal_index=new VECTOR_ND<int>(n);else if(diagonal_index->n!=n)diagonal_index->Resize(n);
    for(int i=1;i<=n;i++){
        SPARSE_VECTOR_ND<T>& row=*A(i);
        for(int j=1;j<=row.number_of_active_indices;j++) if(row.indices[j]==i){(*diagonal_index)(i)=j;continue;}}}

    void Multiply(const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const
    {for(int i=1;i<=n;i++) result(i)=A(i)->Dot_Product(x);}

    void Negate()
    {for(int i=1;i<=n;i++) A(i)->Negate();}

    SPARSE_MATRIX_NXN<T>& operator*=(const T a)
    {for(int i=1;i<=n;i++)(*A(i))*=a;return *this;}

    SPARSE_MATRIX_NXN<T>& operator+=(const T a)
    {for(int i=1;i<=n;i++)Set_Element(i,i,(*this)(i,i)+a);return *this;}

    bool Symmetric(const T tolerance=1e-7)
    {for(int i=1;i<=n;i++){
        SPARSE_VECTOR_ND<T>& row=*A(i);
        for(int j=1;j<=row.number_of_active_indices;j++) if(abs((*this)(row.indices[j],i)-row.x[j])>tolerance) return false;}
    return true;}

    bool Positive_Diagonal_And_Nonnegative_Row_Sum(T tolerance=1e-7)
    {bool return_value=true;
    for(int i=1;i<=n;i++){
        SPARSE_VECTOR_ND<T>& row=*A(i);
        if(row(i)<=0){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            LOG::cout<<"Diagonal Entry "<<i<<" Contains Non Positive Element: "<<row(i)<<std::endl<<row;
#endif
            return false;}
        T sum=0;for(int j=1;j<=row.number_of_active_indices;j++) sum+=row.x[j];
        if(sum<-tolerance){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            LOG::cout<<"Sum Of Row "<<i<<" Is Negative: "<<sum<<std::endl<<row;
#endif
            return_value=false;}}
    return return_value;}

    void Transpose(SPARSE_MATRIX_NXN<T>& A_transpose)
    {assert(A_transpose.n==n);
    for(int i=1;i<=n;i++){SPARSE_VECTOR_ND<T>* row=A(i);for(int j=1;j<=row->number_of_active_indices;j++) A_transpose.Set_Element(row->indices[j],i,row->x[j]);}}

    bool Is_Transpose(SPARSE_MATRIX_NXN<T>& A_transpose,T tolerance)
    {assert(A_transpose.n==n);
    for(int i=1;i<=n;i++) for(int j=1;j<=n;j++) if(abs(A_transpose(i,j)-(*this)(j,i))>tolerance) return false;
    return true;}

    void Solve_Forward_Substitution(const VECTOR_ND<T>& b,VECTOR_ND<T>& x,const bool diagonal_is_identity=false,const bool diagonal_is_inverted=false)
    {if(diagonal_is_identity) for(int i=1;i<=n;i++){
        SPARSE_VECTOR_ND<T>& row=*A(i);
        T sum=0;for(int j=1;j<=(*diagonal_index)(i)-1;j++){sum+=row.x[j]*x(row.indices[j]);}
        x(i)=b(i)-sum;}
    else if(!diagonal_is_inverted) for(int i=1;i<=n;i++){
        SPARSE_VECTOR_ND<T>& row=*A(i);
        T sum=0;for(int j=1;j<=(*diagonal_index)(i)-1;j++){sum+=row.x[j]*x(row.indices[j]);}
        x(i)=(b(i)-sum)/row.x[(*diagonal_index)(i)];}
    else for(int i=1;i<=n;i++){
        SPARSE_VECTOR_ND<T>& row=*A(i);
        T sum=0;for(int j=1;j<=(*diagonal_index)(i)-1;j++){sum+=row.x[j]*x(row.indices[j]);}
        x(i)=(b(i)-sum)*row.x[(*diagonal_index)(i)];}}

    void Solve_Backward_Substitution(const VECTOR_ND<T>& b,VECTOR_ND<T>& x,const bool diagonal_is_identity=false,const bool diagonal_is_inverted=false)
    {if(diagonal_is_identity) for(int i=n;i>=1;i--){
        SPARSE_VECTOR_ND<T>& row=*A(i);
        T sum=0;for(int j=(*diagonal_index)(i)+1;j<=row.number_of_active_indices;j++){sum+=row.x[j]*x(row.indices[j]);}
        x(i)=b(i)-sum;}
    else if(!diagonal_is_inverted) for(int i=n;i>=1;i--){
        SPARSE_VECTOR_ND<T>& row=*A(i);
        T sum=0;for(int j=(*diagonal_index)(i)+1;j<=row.number_of_active_indices;j++){sum+=row.x[j]*x(row.indices[j]);}
        x(i)=(b(i)-sum)/row.x[(*diagonal_index)(i)];}
    else for(int i=n;i>=1;i--){
        SPARSE_VECTOR_ND<T>& row=*A(i);
        T sum=0;for(int j=(*diagonal_index)(i)+1;j<=row.number_of_active_indices;j++){sum+=row.x[j]*x(row.indices[j]);}
        x(i)=(b(i)-sum)*row.x[(*diagonal_index)(i)];}}

    // actually an LU saving square roots, with an inverted diagonal saving divides
    void Construct_Incomplete_Cholesky_Factorization(const bool modified_version=true,const T modified_coefficient=.97,const T zero_tolerance=1e-8,const T zero_replacement=1e-8)
    {delete C;C=new SPARSE_MATRIX_NXN<T>(*this);C->Initialize_Diagonal_Index();
    for(int i=1;i<=n;i++){ // for each row
        SPARSE_VECTOR_ND<T>& row=*C->A(i);int row_diagonal_index=(*C->diagonal_index)(i);T sum=0;
        for(int k_bar=1;k_bar<row_diagonal_index;k_bar++){ // for all the columns before the diagonal element
            int k=row.indices[k_bar];SPARSE_VECTOR_ND<T>& row2=*C->A(k);int row2_diagonal_index=(*C->diagonal_index)(k);
            row.x[k_bar]*=row2.x[row2_diagonal_index]; // divide by the diagonal element (which has already been inverted)
            int j_bar=k_bar+1; // start with the next element in the row, when subtracting the dot product
            for(int i_bar=row2_diagonal_index+1;i_bar<=row2.number_of_active_indices;i_bar++){ // run through the rest of the elements in the row2
                int i=row2.indices[i_bar];T dot_product_term=row.x[k_bar]*row2.x[i_bar];
                while(j_bar<row.number_of_active_indices && row.indices[j_bar]<i) j_bar++; // gets j_bar such that j_bar>=i
                if(row.indices[j_bar]==i) row.x[j_bar]-=dot_product_term;else if(modified_version) sum+=dot_product_term;}}
        T denominator=row.x[row_diagonal_index]-modified_coefficient*sum;
        if(i==n && denominator<=zero_tolerance) row.x[row_diagonal_index]=1/zero_replacement; // ensure last diagonal element is not zero
        else row.x[row_diagonal_index]=1/denominator;}} // finally, store the diagonal element in inverted form

    void Gauss_Seidel_Single_Iteration(VECTOR_ND<T>& x,const VECTOR_ND<T>& b)
    {assert(x.n==b.n && x.n==n);
    for(int i=1;i<=n;i++){
        T rho=0;T diag_entry=0;
        for(int j=1;j<=A(i)->number_of_active_indices;j++){
            int index=A(i)->indices[j];
            if(index==i){diag_entry=A(i)->x[j];continue;}
            rho+=A(i)->x[j]*x(index);}
        x(i)=(b(i)-rho)/diag_entry;}}

    void Gauss_Seidel_Solve(VECTOR_ND<T>& x,const VECTOR_ND<T>& b,const T tolerance=1e-12,const int max_iterations=1000000)
    {assert(x.n==b.n && x.n==n);
    VECTOR_ND<T> last_x(x);
    for(int k=1;k<=max_iterations;k++){
        Gauss_Seidel_Single_Iteration(x,b);
        T residual=0;for(int j=1;j<=n;j++){residual+=sqr(last_x(j)-x(j));last_x(j)=x(j);}if(residual<tolerance) return;}}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    void Write_Row_Lengths()
    {for(int i=1;i<=n;i++)LOG::cout<<A(i)->number_of_active_indices<<" ";LOG::cout<<std::endl;}
#endif

//#####################################################################
};
template<class T> inline std::ostream& operator<<(std::ostream& output_stream,const SPARSE_MATRIX_NXN<T>& A)
{for(int i=1;i<=A.n;i++){for(int j=1;j<=A.n;j++)output_stream<<A(i,j)<<" ";output_stream<<std::endl;}return output_stream;}
}
#endif
