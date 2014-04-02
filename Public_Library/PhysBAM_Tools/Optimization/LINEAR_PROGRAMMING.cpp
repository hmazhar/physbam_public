//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/givens_rotate.h>
#include <PhysBAM_Tools/Optimization/LINEAR_PROGRAMMING.h>
using namespace PhysBAM;
//####################################################################################
// Function Move_Column_From_B_To_N_And_Shift_Down
//####################################################################################
// puts column of B in N, shifts remaining columns of B left, and puts column from S in B
template<class T> void LINEAR_PROGRAMMING<T>::
Move_Column_From_B_To_N_And_Shift_Down(MATRIX_MXN<T>& B,VECTOR_ND<int>& permute_B,const int b_column,MATRIX_MXN<T>& N,VECTOR_ND<int>& permute_N,const int n_column)
{
    assert(B.m==B.n);
    for(int i=1;i<=B.n;i++)N(i,n_column)=B(i,b_column);permute_N(n_column)=permute_B(b_column);
    for(int j=b_column;j<B.n;j++){for(int i=1;i<=B.n;i++)B(i,j)=B(i,j+1);permute_B(j)=permute_B(j+1);}
}
//####################################################################################
// Function Update_Upper_Triangular_Matrix_After_Column_Shift
//####################################################################################
// columns starting at given index have been shifted left so they have one value below diagonal
// but we can make use of the fact that those values used to be on the diagonal (must be good values)
// also updates N and b with the transformation
template<class T> void LINEAR_PROGRAMMING<T>::
Update_Upper_Triangular_Matrix_After_Column_Shift(MATRIX_MXN<T>& A,MATRIX_MXN<T>& N,VECTOR_ND<T>& b,const int column,const T tolerance,const bool check_last_column)
{
    assert(A.m==A.n);
    for(int i=column;i<A.n;i++){
        assert(abs(A(i+1,i))>tolerance);
        VECTOR<T,2> v(A(i,i),A(i+1,i));A(i,i)=v.Magnitude();A(i+1,i)=0;v.Normalize();
        for(int j=i+1;j<=A.n;j++) givens_rotate(A(i,j),A(i+1,j),v.x,v.y);
        for(int j=1;j<=N.n;j++) givens_rotate(N(i,j),N(i+1,j),v.x,v.y);
        givens_rotate(b(i),b(i+1),v.x,v.y);}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(check_last_column && abs(A(A.n,A.n))<tolerance){
        LOG::cout << "Found near-zero on diagonal during attempt to update upper triangular matrix!" << std::endl;
        LOG::cout << A << std::endl;
        /*PHYSBAM_FATAL_ERROR()*/}
#endif
}
//####################################################################################
// Function Find_Feasible_Solution
//####################################################################################
template<class T> void LINEAR_PROGRAMMING<T>::
Find_Feasible_Solution(MATRIX_MXN<T>& B,MATRIX_MXN<T>& N,VECTOR_ND<T>& x_B,VECTOR_ND<T>& b,VECTOR_ND<T>& x_N,VECTOR_ND<int>& permute_B,VECTOR_ND<int>& permute_N,VECTOR_ND<T>& x_unpermuted,
    const ARRAY<PAIR<bool,T> >& x_min,const ARRAY<PAIR<bool,T> >& x_max,const T tolerance,const T step_tolerance,bool verbose)
{
    assert(B.m==B.n);
    ARRAY<bool> x_in_feasible_region(B.n+N.n);
    for(int i=1;i<=x_N.n;i++){ // TODO: starting x_N at min values (or max if no min)
        if(x_min(permute_N(i)).x) x_N(i)=x_min(permute_N(i)).y;
        else{assert(x_max(permute_N(i)).x);x_N(i)=x_max(permute_N(i)).y;} // must have a max if no min
        x_in_feasible_region(permute_N(i))=true;}

    VECTOR_ND<T> c_B(B.n);
    int iteration=1;
    for(;;iteration++){
        // solve
        x_B=B.Upper_Triangular_Solve(b-N*x_N);

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(verbose){
            LOG::cout<<"************* ITERATION "<<iteration<< ":\nB:\n"<<B<<"\nN:\n"<<N<<"\nb:\n"<<b<<"\nx_N:\n"<<x_N<<std::endl;
            LOG::cout<<"permute_B:\n"<<permute_B<<"\npermute_N:\n"<<permute_N<<"\nx_B:\n"<<x_B<<std::endl;
            MATRIX_MXN<T> inverse;B.LU_Inverse(inverse);LOG::cout<<"B^-1 * N:\n"<<inverse*N<<"\nB^-1 * b:\n"<<inverse*b<<std::endl;}
#endif

        bool have_unsatisfied_constraint=false;
        // TODO: maybe use tolerance to consider anything in [x_min-tol,x_max+tol] as being in feasible region
        for(int i=1;i<=B.n;i++){
            if(!x_in_feasible_region(permute_B(i)) && x_min(permute_B(i)).x && x_B(i)<x_min(permute_B(i)).y){c_B(i)=-1;have_unsatisfied_constraint=true;}
            else if(!x_in_feasible_region(permute_B(i)) && x_max(permute_B(i)).x && x_B(i)>x_max(permute_B(i)).y){c_B(i)=+1;have_unsatisfied_constraint=true;}
            else{ // clamp ones we know should be in feasible region
                x_in_feasible_region(permute_B(i))=true;c_B(i)=0;
                if(x_min(permute_B(i)).x) x_B(i)=max(x_B(i),x_min(permute_B(i)).y);
                if(x_max(permute_B(i)).x) x_B(i)=min(x_B(i),x_max(permute_B(i)).y);}}

        if(!have_unsatisfied_constraint){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(verbose)LOG::cout << "All constraints satisfied!" << std::endl;
#endif
            break;} // done

        VECTOR_ND<T> pi=B.Transpose_Lower_Triangular_Solve(c_B),sigma(-N.Transpose_Times(pi));

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(verbose) LOG::cout<<"c_B =\n"<<c_B<<"pi =\n"<<pi<<"sigma =\n"<<sigma<<std::endl;
#endif

        int index_to_release=0;T sigma_to_release=0;
        for(int i=1;i<=N.n;i++) if((sigma(i)<0 && x_min(permute_N(i)).x && x_N(i)==x_min(permute_N(i)).y) ||
                                   (sigma(i)>0 && x_max(permute_N(i)).x && x_N(i)==x_max(permute_N(i)).y))
            if(abs(sigma(i))>sigma_to_release){sigma_to_release=abs(sigma(i));index_to_release=i;}
      
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        T sigma_tolerance=(T)1e-10;
        if(sigma_to_release<sigma_tolerance || !index_to_release){
            LOG::cout << "No feasible point found for given constraints!" << std::endl;
            break;} // TODO: how handle this case?

        if(verbose) LOG::cout << "Releasing index " << index_to_release << " (muscle " << permute_N(index_to_release) << ")" << std::endl;
#endif

        VECTOR_ND<T> p_N(N.n);
        if(x_min(permute_N(index_to_release)).x && x_N(index_to_release)==x_min(permute_N(index_to_release)).y) p_N(index_to_release)=1;
        else if(x_max(permute_N(index_to_release)).x && x_N(index_to_release)==x_max(permute_N(index_to_release)).y) p_N(index_to_release)=-1;
        VECTOR_ND<T> p_B(-B.Upper_Triangular_Solve(N*p_N));

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(verbose) LOG::cout<<"Got p_B\n"<<p_B<<"\nGot p_N\n"<<p_N<<std::endl;
#endif

        T alpha=FLT_MAX,limiting_value=0;int limiting_index=0;
        // clamp based on other constraint for variable we are releasing (limiting_index=-1 used to indicate this case)
        if(x_min(permute_N(index_to_release)).x && x_max(permute_N(index_to_release)).x){
            limiting_value=(p_N(index_to_release)>0)?x_max(permute_N(index_to_release)).y:x_min(permute_N(index_to_release)).y;
            alpha=p_N(index_to_release)*(limiting_value-x_N(index_to_release));
            limiting_index=-1;}
        // check if we hit any other constraint
        // NOTE: if p_B(i) is close to zero we don't clamp against that variable because if we did and ended up trying to swap
        // columns from N and B then the new B would not be full rank and trying to update it using LU would give us a zero on the diagonal
        for(int i=1;i<=B.n;i++) if(!x_in_feasible_region(permute_B(i)) || abs(p_B(i)*alpha)>=step_tolerance){
            if(p_B(i)>0 && x_max(permute_B(i)).x && (x_in_feasible_region(permute_B(i)) || (x_max(permute_B(i)).y-x_B(i))>0) && (x_max(permute_B(i)).y-x_B(i))/p_B(i)<alpha){
                alpha=max((T)0,(x_max(permute_B(i)).y-x_B(i))/p_B(i));limiting_index=i;limiting_value=x_max(permute_B(i)).y;}
            else if(p_B(i)<0 && x_min(permute_B(i)).x && (x_in_feasible_region(permute_B(i)) || (x_min(permute_B(i)).y-x_B(i))<0) && (x_min(permute_B(i)).y-x_B(i))/p_B(i)<alpha){
                alpha=max((T)0,(x_min(permute_B(i)).y-x_B(i))/p_B(i));limiting_index=i;limiting_value=x_min(permute_B(i)).y;}}

        assert(alpha>=0);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(verbose) LOG::cout<<"alpha="<<alpha<<", limiting_index="<<limiting_index<<", limiting_value="<<limiting_value<<"\nx_in_feasible_region = "<<x_in_feasible_region<<std::endl;
#endif

        // robustly detect when a variable enters feasible region as a result of moving by alpha
        for(int i=1;i<=B.n;i++) if(!x_in_feasible_region(permute_B(i))){
            if((p_B(i)>0 && x_min(permute_B(i)).x && x_B(i)<=x_min(permute_B(i)).y && (x_min(permute_B(i)).y-x_B(i))/p_B(i)<=alpha) ||
               (p_B(i)<0 && x_max(permute_B(i)).x && x_B(i)>=x_max(permute_B(i)).y && (x_max(permute_B(i)).y-x_B(i))/p_B(i)<=alpha))
                x_in_feasible_region(permute_B(i))=true;}

        if(!limiting_index){
#if 0
            alpha=0;limiting_index=0; // see how far we need to go to satisfy as much as possible
            for(int i=1;i<=B.n;i++)
                if(p_B(i)>0 && x_min(permute_B(i)).x && (x_min(permute_B(i)).y-x_B(i))/p_B(i)>alpha){alpha=(x_min(permute_B(i)).y-x_B(i))/p_B(i);limiting_index=i;}
                else if(p_B(i)<0 && x_max(permute_B(i)).x && (x_max(permute_B(i)).y-x_B(i))/p_B(i)>alpha){alpha=(x_max(permute_B(i)).y-x_B(i))/p_B(i);limiting_index=i;}
            x_B+=alpha*p_B;
#endif
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            LOG::cout << "No limiting index, picking alpha="<<alpha<<" (reaching index"<<limiting_index<<")"<<std::endl;
#endif
            break;} // TODO: how handle this case?

        if(limiting_index==-1){ // p_N(i) limits alpha, switch x_N to be opposite end's constraint
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(verbose) LOG::cout << "Active constraint " << index_to_release << " (muscle " << permute_N(index_to_release) << ") switched to bound " << limiting_value << std::endl;
#endif
            x_N(index_to_release)=limiting_value;}
        else{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(verbose){
                LOG::cout << "limiting_index="<<limiting_index<<" (muscle " << permute_B(limiting_index) << ") getting constrained to " << limiting_value 
                          << ", index_to_release="<<index_to_release<<" (muscle " << permute_N(index_to_release) << ") released from " << x_N(index_to_release) << ")" << std::endl;
                LOG::cout << "BEFORE COLUMN SWITCH ("<<limiting_index<<","<<index_to_release<<")\nB:\n"<<B<<"\nN:\n"<<N<<std::endl;}
#endif

            // put B column "limiting_index" in N column "index_to_release", and put N column "index_to_release" in
            // the last column of B (after shifting the other columns left)
            int permute_n=permute_N(index_to_release);VECTOR_ND<T> column_from_N(B.n);N.Get_Column(index_to_release,column_from_N);
            Move_Column_From_B_To_N_And_Shift_Down(B,permute_B,limiting_index,N,permute_N,index_to_release);
            B.Set_Column(B.n,column_from_N);permute_B(B.n)=permute_n;

            // update x_N
            x_N(index_to_release)=limiting_value;
            x_in_feasible_region(permute_N(index_to_release))=true;

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(verbose) LOG::cout << "AFTER COLUMN SWITCH\nB:\n"<<B<<"\nN:\n"<<N<<std::endl;
#endif

            Update_Upper_Triangular_Matrix_After_Column_Shift(B,N,b,limiting_index,tolerance); // make upper triangular again
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(verbose) LOG::cout << "AFTER MAKING UPPER TRIANGULAR\nB:\n"<<B<<"\nN:\n"<<N<<std::endl;
#endif
        }
    }

    if(x_in_feasible_region.Number_False()){ // have some values not in feasible region 
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout << "PHYSBAM_WARNING: clamping values into feasible region after unsuccessful LP solve" << std::endl;
#endif
        for(int i=1;i<=x_B.n;i++) if(!x_in_feasible_region(permute_B(i))){
            if(x_min(permute_B(i)).x) x_B(i)=max(x_B(i),x_min(permute_B(i)).y);
            if(x_max(permute_B(i)).x) x_B(i)=min(x_B(i),x_max(permute_B(i)).y);}}

    VECTOR_ND<int> permute(B.n+N.n);permute.Set_Subvector(1,permute_B);permute.Set_Subvector(B.n+1,permute_N);
    VECTOR_ND<T> x(B.n+N.n);x.Set_Subvector(1,x_B);x.Set_Subvector(B.n+1,x_N);
    x_unpermuted=x.Unpermute(permute);

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    LOG::cout << "LP finished in " << iteration << " iterations" << std::endl;
    if(verbose) LOG::cout<<"LP Result x_B:\n"<<x_B<<"\nLP Result permute:\n"<<permute<<std::endl;
#endif
}
//####################################################################################
template class LINEAR_PROGRAMMING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_PROGRAMMING<double>;
#endif
