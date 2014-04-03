//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUADRATIC_PROGRAMMING
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/givens_rotate.h>
#include <PhysBAM_Tools/Optimization/QUADRATIC_PROGRAMMING.h>
using namespace PhysBAM;
//####################################################################################
// Matrix utilities
//####################################################################################
// puts column of B in N, shifts remaining columns of B left, and puts column from S in B
template<class T> void QUADRATIC_PROGRAMMING<T>::
Move_Column_From_B_To_N_And_Shift_Down(MATRIX_MXN<T>& B,VECTOR_ND<int>& permute_B,const int b_column,MATRIX_MXN<T>& N,VECTOR_ND<int>& permute_N,const int n_column)
{
    assert(B.m==B.n);
    for(int i=1;i<=B.n;i++)N(i,n_column)=B(i,b_column);permute_N(n_column)=permute_B(b_column);
    for(int j=b_column;j<B.n;j++){for(int i=1;i<=B.n;i++)B(i,j)=B(i,j+1);permute_B(j)=permute_B(j+1);}
}
template<class T> void QUADRATIC_PROGRAMMING<T>::
Remove_Column(const int column,MATRIX_MXN<T>& A,VECTOR_ND<int>& permute,VECTOR_ND<T>* x)
{
    if(column!=A.n){
        for(int i=1;i<=A.m;i++) A(i,column)=A(i,A.n);
        permute(column)=permute(A.n);
        if(x) (*x)(column)=(*x)(A.n);}
    A.Resize(A.m,A.n-1);permute.Resize(A.n);if(x)x->Resize(A.n);
}
//####################################################################################
// Update_Upper_Triangular_Matrix_After_Column_Shift
//####################################################################################
// columns starting at given index have been shifted left so they have one value below diagonal
// but we can make use of the fact that those values used to be on the diagonal (must be good values)
// also updates N and b with the transformation
template<class T> void QUADRATIC_PROGRAMMING<T>::
Update_Upper_Triangular_Matrix_After_Column_Shift(MATRIX_MXN<T>& A,MATRIX_MXN<T>& S,MATRIX_MXN<T>& N,VECTOR_ND<T>& b,const int column,const T tolerance,const bool check_last_column)
{
    assert(A.m==A.n);
    for(int i=column;i<A.n;i++){
        assert(abs(A(i+1,i))>tolerance);
        VECTOR<T,2> v(A(i,i),A(i+1,i));A(i,i)=v.Magnitude();A(i+1,i)=0;v.Normalize();
        for(int j=i+1;j<=A.n;j++) givens_rotate(A(i,j),A(i+1,j),v.x,v.y);
        for(int j=1;j<=S.n;j++) givens_rotate(S(i,j),S(i+1,j),v.x,v.y);
        for(int j=1;j<=N.n;j++) givens_rotate(N(i,j),N(i+1,j),v.x,v.y);
        givens_rotate(b(i),b(i+1),v.x,v.y);}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(check_last_column && abs(A(A.n,A.n))<tolerance){
        LOG::cout << "Found near-zero on diagonal during attempt to update upper triangular matrix!" << std::endl;
        LOG::cout << A << std::endl;
        /*PHYSBAM_FATAL_ERROR();*/}
#endif
}
//####################################################################################
// Find_Optimal_Solution
//####################################################################################
// assumes given variables are in feasible region
template<class T> void QUADRATIC_PROGRAMMING<T>::
Find_Optimal_Solution(MATRIX_MXN<T>& B,MATRIX_MXN<T>& S,MATRIX_MXN<T>& N,VECTOR_ND<T>& x_B,VECTOR_ND<T>& x_S,VECTOR_ND<T>& b,VECTOR_ND<T>& b_N,VECTOR_ND<int>& permute_B,
    VECTOR_ND<int>& permute_S,VECTOR_ND<int>& permute_N,const MATRIX_MXN<T>& D_unpermuted,const MATRIX_MXN<T>& epsilon_hat_unpermuted,const VECTOR_ND<T>& f_hat,
    VECTOR_ND<T>& x_unpermuted,const ARRAY<PAIR<bool,T> >& x_min,const ARRAY<PAIR<bool,T> >& x_max,const T tolerance,const T step_tolerance,const bool debug_optimization)
{
    assert(B.m==B.n);assert(D_unpermuted.m==D_unpermuted.n);
    MATRIX_MXN<T> D(D_unpermuted.n);MATRIX_MXN<T> epsilon_hat(epsilon_hat_unpermuted.m);
    VECTOR_ND<int> permute(B.n+S.n+N.n);
    VECTOR_ND<T> x(B.n+S.n+N.n);

    int iteration=1;
    for(;;iteration++){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(debug_optimization){
            LOG::cout<<"************* ITERATION "<<iteration<< ":\nB:\n"<<B<<"\nS:\n"<<S<<"\nN:\n"<<N<<"\nb:\n"<<b<<"\nb_N:\n"<<b_N<<std::endl;
            LOG::cout<<"permute_B:\n"<<permute_B<<"\npermute_S:\n"<<permute_S<<"\npermute_N:\n"<<permute_N<<"\n\nx_B:\n"<<x_B<<"\nx_S:\n"<<x_S<<std::endl;}
#endif

        // update permute and full x vector
        permute.Set_Subvector(1,permute_B);permute.Set_Subvector(B.n+1,permute_S);permute.Set_Subvector(B.n+S.n+1,permute_N);
        D=D_unpermuted.Permute_Columns(permute);epsilon_hat=epsilon_hat_unpermuted.Permute_Columns(permute);
        x.Set_Subvector(1,x_B);x.Set_Subvector(B.n+1,x_S);x.Set_Subvector(B.n+S.n+1,b_N); // note that b_N==x_N

        // find search direction
        VECTOR_ND<T> p_S,p_B;
        T alpha=1,limiting_value=0;int limiting_index=0; // negative limiting_index will indicate it comes from the "S" set

        // if S matrix is empty (e.g. in first iteration) then the null space is empty and there are no feasible directions.
        // so we skip the step below and go to lagrange multiplier check
        if(x_S.n){
            MATRIX_MXN<T> negative_B_inverse_S=-B.Upper_Triangular_Solve(S);
            MATRIX_MXN<T> Z(B.n+S.n+N.n,S.n),D_hat_times_Z(x.n+f_hat.n,S.n);
            Z.Add_To_Submatrix(1,1,negative_B_inverse_S);for(int i=1;i<=S.n;i++)Z(B.n+i,i)=1;
            D_hat_times_Z.Add_To_Submatrix(1,1,D*Z);D_hat_times_Z.Add_To_Submatrix(x.n+1,1,epsilon_hat*Z);

            VECTOR_ND<T> rhs(x.n+f_hat.n);rhs.Set_Subvector(1,D*(-x));rhs.Set_Subvector(x.n+1,f_hat-epsilon_hat*x);
            VECTOR_ND<T> p_S=D_hat_times_Z.Normal_Equations_Solve(rhs);
            VECTOR_ND<T> p_B=negative_B_inverse_S*p_S;
           
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(debug_optimization) LOG::cout<<"Got p_B\n"<<p_B<<"\nGot p_S\n"<<p_S<<std::endl;
#endif

            // check ones in B
            for(int i=1;i<=B.n;i++) if(abs(p_B(i)*alpha)>=step_tolerance){
                if(p_B(i)>0 && x_max(permute_B(i)).x && (x_max(permute_B(i)).y-x_B(i))/p_B(i)<alpha){
                    alpha=max((T)0,(x_max(permute_B(i)).y-x_B(i))/p_B(i));limiting_index=i;limiting_value=x_max(permute_B(i)).y;}
                else if(p_B(i)<0 && x_min(permute_B(i)).x && (x_min(permute_B(i)).y-x_B(i))/p_B(i)<alpha){
                    alpha=max((T)0,(x_min(permute_B(i)).y-x_B(i))/p_B(i));limiting_index=i;limiting_value=x_min(permute_B(i)).y;}}
            // check ones in S
            for(int i=1;i<=S.n;i++) if(abs(p_S(i)*alpha)>=step_tolerance){
                if(p_S(i)>0 && x_max(permute_S(i)).x && (x_max(permute_S(i)).y-x_S(i))/p_S(i)<alpha){
                    alpha=max((T)0,(x_max(permute_S(i)).y-x_S(i))/p_S(i));limiting_index=-i;limiting_value=x_max(permute_S(i)).y;}
                else if(p_S(i)<0 && x_min(permute_S(i)).x && (x_min(permute_S(i)).y-x_S(i))/p_S(i)<alpha){
                    alpha=max((T)0,(x_min(permute_S(i)).y-x_S(i))/p_S(i));limiting_index=-i;limiting_value=x_min(permute_S(i)).y;}}

            assert(alpha>=0);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(debug_optimization) LOG::cout << "alpha="<<alpha<<", limiting_index="<<limiting_index<<", limiting_value="<<limiting_value<<std::endl;
#endif

            x_B+=alpha*p_B;x_S+=alpha*p_S;

            // reconstruct (permuted) x vector (note that b_N==x_N)
            x.Set_Subvector(1,x_B);x.Set_Subvector(B.n+1,x_S);x.Set_Subvector(B.n+S.n+1,b_N);}

        if(!limiting_index){assert(alpha==1); // check lagrange multipliers
            VECTOR_ND<T> gradient(D.Transpose_Times(D)*x+epsilon_hat.Transpose_Times(epsilon_hat*x-f_hat));

            // decompose gradient
            VECTOR_ND<T> gradient_B(B.n),gradient_N(N.n);gradient.Get_Subvector(1,gradient_B);gradient.Get_Subvector(B.n+S.n+1,gradient_N);
            VECTOR_ND<T> pi=B.Transpose_Lower_Triangular_Solve(gradient_B);
            VECTOR_ND<T> sigma(gradient_N-N.Transpose_Times(pi));

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(debug_optimization) LOG::cout<<"gradient_B =\n"<<gradient_B<<"\ngradient_N =\n"<<gradient_N<<"\npi =\n"<<pi<<"\nsigma =\n"<<sigma<<std::endl;
#endif
            
            int index_to_release=0;T sigma_to_release=0;
            for(int i=1;i<=N.n;i++) if((sigma(i)<0 && x_min(permute_N(i)).x && b_N(i)==x_min(permute_N(i)).y) ||
                                       (sigma(i)>0 && x_max(permute_N(i)).x && b_N(i)==x_max(permute_N(i)).y))
                if(abs(sigma(i))>sigma_to_release){sigma_to_release=abs(sigma(i));index_to_release=i;}
            
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(!index_to_release){if(debug_optimization)LOG::cout << "QP: Nothing to release!" << std::endl;break;} // TODO: how handle this case?

            if(debug_optimization) LOG::cout << "Releasing index " << index_to_release << " (muscle " << permute_N(index_to_release) << ")" << std::endl;
#endif
            
            // move column of N (index_to_release) into S (update S and permute_S)
            S.Resize(S.m,S.n+1);x_S.Resize(S.n);permute_S.Resize(S.n);
            for(int i=1;i<=S.m;i++) S(i,S.n)=N(i,index_to_release);
            x_S(S.n)=b_N(index_to_release);
            permute_S(S.n)=permute_N(index_to_release);

            Remove_Column(index_to_release,N,permute_N,&b_N);}
        else if(limiting_index<0){ // constrain variable from S
            limiting_index*=-1; // make positive

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(debug_optimization) LOG::cout << "Constraining S index " << limiting_index << " (muscle " << permute_S(limiting_index) << ") to " << limiting_value << std::endl;
#endif

            // move column of S (limiting_index) into N (update N, b_N, and permute_N)
            N.Resize(N.m,N.n+1);b_N.Resize(N.n);permute_N.Resize(N.n);
            for(int i=1;i<=N.m;i++) N(i,N.n)=S(i,limiting_index);
            b_N(N.n)=limiting_value;
            permute_N(N.n)=permute_S(limiting_index);

            Remove_Column(limiting_index,S,permute_S,&x_S);}
        else{ // constrain variable from B
            assert(S.n); // S should not be empty if we got to this case (we'd be doing lagrange multipliers instead)

            if(debug_optimization){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
                LOG::cout << "Constraining B index " << limiting_index << " (muscle " << permute_B(limiting_index) << ") to " << limiting_value << std::endl;
                LOG::cout << "BEFORE COLUMN SWITCH\nB:\n" << B << "\nS:\n" << S << "\nN:\n" << N << std::endl;
#endif
            }

            // move column from B to N
            N.Resize(N.m,N.n+1);b_N.Resize(N.n);permute_N.Resize(N.n);
            Move_Column_From_B_To_N_And_Shift_Down(B,permute_B,limiting_index,N,permute_N,N.n);
            b_N(N.n)=limiting_value;for(int i=limiting_index;i<B.n;i++)x_B(i)=x_B(i+1);

            // re-factor B (ignoring last column, which we will attempt to fill from S next)
            Update_Upper_Triangular_Matrix_After_Column_Shift(B,S,N,b,limiting_index,tolerance,false); // make upper triangular again
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(debug_optimization) LOG::cout << "AFTER COLUMN B->N AND REFACTORING B\nB:\n" << B << "\nS:\n" << S << "\nN:\n" << N << std::endl;
#endif

            // now choose the best column from S to stick in last column of B
            int s_column=0;T largest_s_diagonal=0;
            for(int j=1;j<=S.n;j++)if(abs(S(S.m,j))>largest_s_diagonal){s_column=j;largest_s_diagonal=abs(S(S.m,j));}
            if(!s_column){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
                LOG::cout<<"Bottom row of B,S is all zero!"<<std::endl;
#endif
                break;}

            // copy column from S to B and remove from S
            for(int i=1;i<=B.n;i++)B(i,B.n)=S(i,s_column);permute_B(B.n)=permute_S(s_column);x_B(B.n)=x_S(s_column);
            Remove_Column(s_column,S,permute_S,&x_S);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(debug_optimization) LOG::cout << "AFTER MOVING S COLUMN ("<<s_column<<") TO B\nB:\n" << B << "\nS:\n" << S << "\nN:\n" << N << std::endl;
#endif
        }
    }

    x.Set_Subvector(1,x_B);x.Set_Subvector(B.n+1,x_S);x.Set_Subvector(B.n+S.n+1,b_N); // note that b_N==x_N
    permute.Set_Subvector(1,permute_B);permute.Set_Subvector(B.n+1,permute_S);permute.Set_Subvector(B.n+S.n+1,permute_N);
    x_unpermuted=x.Unpermute(permute);

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(debug_optimization) LOG::cout << "Result x_B:\n"<<x_B<<"\nResult x_S:\n"<<x_S<<"\nResult permute:\n"<<permute<<"\nResult x:\n"<<x<<std::endl;

    LOG::cout << "QP finished in " << iteration << " iterations" << std::endl;
#endif
}
//####################################################################################
template class QUADRATIC_PROGRAMMING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class QUADRATIC_PROGRAMMING<double>;
#endif
