//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/INT_ITERATOR_THREADED.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_THREADED.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>
using namespace PhysBAM;
//#####################################################################
// Function Parallel_Solve
//#####################################################################
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance)
{
    Init_Barriers();
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    assert(tid==domain_index(interior_domain.max_corner));
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);
    const ARRAY<INTERVAL<int> >& ghost_indices=all_ghost_indices(tid);

    int global_n=A.n,interior_n=interior_indices.Size()+1;
    T global_tolerance=tolerance;
    int desired_iterations=global_n;if(enforce_compatibility) desired_iterations--;if(maximum_iterations) desired_iterations=min(desired_iterations,maximum_iterations);

    VECTOR_ND<T> z_interior(interior_n,false);
#ifdef USE_PTHREADS
    pthread_mutex_lock(&sum_lock);
#endif
    temp.Resize(global_n);p.Resize(global_n);
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&sum_lock);
#endif

    // build interior views of x,b,p,z,temp
    VECTOR_ND<T> x_interior,b_interior,p_interior,temp_interior;
    x_interior.Set_Subvector_View(x,interior_indices);
    b_interior.Set_Subvector_View(b,interior_indices);
    p_interior.Set_Subvector_View(p,interior_indices);
    temp_interior.Set_Subvector_View(temp,interior_indices);

    // adjust x for the null space
    if(enforce_compatibility && remove_null_space_solution_component) x_interior-=(T)(Global_Sum(x_interior.Sum_Double_Precision(),tid)/global_n);

    // find initial residual, r=b-Ax - reusing b for the residual
#ifdef USE_PTHREADS
    pthread_barrier_wait(&barr);
#endif
    A.Times(interior_indices,ghost_indices,x,temp);
    b_interior-=temp_interior;
    if(enforce_compatibility) b_interior-=(T)(Global_Sum(b_interior.Sum_Double_Precision(),tid)/global_n);
    if(Global_Max(b_interior.Max_Abs())<=global_tolerance){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(show_results) LOG::cout<<"NO ITERATIONS NEEDED"<<std::endl;
#endif
        return;}

    SPARSE_MATRIX_FLAT_NXN<T>* C=0;
    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
    if(incomplete_cholesky){
        C=A.Create_Submatrix(interior_indices);
        C->In_Place_Incomplete_Cholesky_Factorization(modified_incomplete_cholesky,modified_incomplete_cholesky_coefficient,
            preconditioner_zero_tolerance,preconditioner_zero_replacement);}

    double rho=0,rho_old=0;
    for(int iteration=1;;iteration++){
        if(incomplete_cholesky){
            // solve Mz=r
            C->Solve_Forward_Substitution(b_interior,temp_interior,true); // diagonal should be treated as the identity
            C->Solve_Backward_Substitution(temp_interior,z_interior,false,true);} // diagonal is inverted to save on divides
        else z_interior=b_interior; // set z=r when there is no preconditioner

        // for Neumann boundary conditions only, make sure z sums to zero
        if(enforce_compatibility) z_interior-=(T)(Global_Sum(z_interior.Sum_Double_Precision(),tid)/global_n);

        // update search direction
        rho_old=rho;rho=Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(z_interior,b_interior),tid);
        T beta=0;if(iteration==1) p_interior=z_interior;else{beta=(T)(rho/rho_old);for(int i=1;i<=interior_n;i++) p_interior(i)=z_interior(i)+beta*p_interior(i);} // when iteration=1, beta=0

        // update solution and residual
#ifdef USE_PTHREADS
        pthread_barrier_wait(&barr);
#endif
        A.Times(interior_indices,ghost_indices,p,temp);
        T alpha=(T)(rho/Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(p_interior,temp_interior),tid));
        for(int i=1;i<=interior_n;i++){x_interior(i)+=alpha*p_interior(i);b_interior(i)-=alpha*temp_interior(i);}

        // remove null space component of b before computing residual norm because we might have converged up to the null space but have some null space component left due to roundoff
        if(enforce_compatibility) b_interior-=(T)(Global_Sum(b_interior.Sum_Double_Precision(),tid)/global_n);

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        T residual=Global_Max(b_interior.Max_Abs());

        // check for convergence
        if(show_residual) LOG::cout<<residual<<std::endl;
        if(residual<=global_tolerance){if(show_results) LOG::cout<<"NUMBER OF ITERATIONS = "<<iteration<<std::endl;break;}
        if(iteration==desired_iterations){if(show_results) LOG::cout<<"DID NOT CONVERGE IN "<<iteration<<" ITERATIONS"<<std::endl;break;}
#endif
    }
    delete C;
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_In_Parts(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance)
{
    int global_n=A.n;
    T global_tolerance=tolerance;
    int desired_iterations=global_n;if(enforce_compatibility) desired_iterations--;if(maximum_iterations) desired_iterations=min(desired_iterations,maximum_iterations);
    VECTOR_ND<T> z(global_n,false);
    temp.Resize(global_n);p.Resize(global_n);
    
    INT_ITERATOR_THREADED_ALPHA<PCG_SPARSE_THREADED<TV> > threaded_iterator(1,global_n,&thread_queue);
    int num_intervals=threaded_iterator.intervals.m;
    ARRAY<T> local_sum(num_intervals);
    
    // adjust x for the null space
    T sum=0;
    if(enforce_compatibility && remove_null_space_solution_component){
        threaded_iterator.template Run<VECTOR_ND<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Sum,x,local_sum);
        for(int i=1;i<=num_intervals;i++) sum+=local_sum(i);}
    
    // find initial residual, r=b-Ax - reusing b for the residual
    threaded_iterator.template Run<SPARSE_MATRIX_FLAT_NXN<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,T>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Part_One,A,x,b,sum/global_n);
    if(enforce_compatibility){sum=0;
        threaded_iterator.template Run<VECTOR_ND<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Sum,b,local_sum);
        for(int i=1;i<=num_intervals;i++) sum+=local_sum(i);
        threaded_iterator.template Run<VECTOR_ND<T>&,T>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Subtract,b,sum/global_n);}

    threaded_iterator.template Run<VECTOR_ND<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Max,b,local_sum);
    T max_val=0;for(int i=1;i<=num_intervals;i++) max_val=max(max_val,local_sum(i));
    if(max_val<=global_tolerance){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(show_results) LOG::cout<<"NO ITERATIONS NEEDED"<<std::endl;
#endif
        return;}

    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
    SPARSE_MATRIX_FLAT_NXN<T>* C=0;
    if(incomplete_cholesky){
        C=new SPARSE_MATRIX_FLAT_NXN<T>(A);
        C->In_Place_Incomplete_Cholesky_Factorization(modified_incomplete_cholesky,modified_incomplete_cholesky_coefficient,
            preconditioner_zero_tolerance,preconditioner_zero_replacement);}

    double rho=0,rho_old=0;
    for(int iteration=1;;iteration++){
       if(incomplete_cholesky){
            // solve Mz=r
            C->Solve_Forward_Substitution(b,temp,true); // diagonal should be treated as the identity
            C->Solve_Backward_Substitution(temp,z,false,true);} // diagonal is inverted to save on divides
        else z=b; // set z=r when there is no preconditioner

        // for Neumann boundary conditions only, make sure z sums to zero
        if(enforce_compatibility){T sum=0;
            threaded_iterator.template Run<VECTOR_ND<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Sum,z,local_sum);
            for(int i=1;i<=num_intervals;i++) sum+=local_sum(i);
            threaded_iterator.template Run<VECTOR_ND<T>&,T>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Subtract,z,sum/global_n);}
        
        // update search direction
        threaded_iterator.template Run<VECTOR_ND<T>&,VECTOR_ND<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Dot,z,b,local_sum);
        rho_old=rho;rho=0;for(int i=1;i<=num_intervals;i++) rho+=local_sum(i);
        threaded_iterator.template Run<VECTOR_ND<T>&,T,T,int>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Part_Two,z,rho,rho_old,iteration);
        
        // update solution and residual
        threaded_iterator.template Run<SPARSE_MATRIX_FLAT_NXN<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Part_Three,A);
        threaded_iterator.template Run<VECTOR_ND<T>&,VECTOR_ND<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Dot,p,temp,local_sum);
        T sum=0;for(int i=1;i<=num_intervals;i++) sum+=local_sum(i);
        threaded_iterator.template Run<VECTOR_ND<T>&,VECTOR_ND<T>&,T>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Part_Four,x,b,rho/sum);

        // remove null space component of b before computing residual norm because we might have converged up to the null space but have some null space component left due to roundoff
        if(enforce_compatibility){T sum=0;
            threaded_iterator.template Run<VECTOR_ND<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Sum,b,local_sum);
            for(int i=1;i<=num_intervals;i++) sum+=local_sum(i);
            threaded_iterator.template Run<VECTOR_ND<T>&,T>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Subtract,b,sum/global_n);}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        threaded_iterator.template Run<VECTOR_ND<T>&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Threaded_Max,b,local_sum);
        T residual=0;for(int i=1;i<=num_intervals;i++) residual=max(residual,local_sum(i));

        // check for convergence
        if(show_residual) LOG::cout<<residual<<std::endl;
        if(residual<=global_tolerance){if(show_results) LOG::cout<<"NUMBER OF ITERATIONS = "<<iteration<<std::endl;break;}
        if(iteration==desired_iterations){if(show_results) LOG::cout<<"DID NOT CONVERGE IN "<<iteration<<" ITERATIONS"<<std::endl;break;}
#endif
    }
    delete C;
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_In_Parts(DOMAIN_ITERATOR_THREADED_ALPHA<PCG_SPARSE_THREADED<TV>,TV>& threaded_iterator,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance)
{
    int global_n=A.n,num_domains=threaded_iterator.domains.m;
    T global_tolerance=tolerance;
    int desired_iterations=global_n;if(enforce_compatibility) desired_iterations--;if(maximum_iterations) desired_iterations=min(desired_iterations,maximum_iterations);
    ARRAY<VECTOR_ND<T> > z_interior(num_domains),x_interior(num_domains),b_interior(num_domains),p_interior(num_domains),temp_interior(num_domains);
    temp.Resize(global_n);p.Resize(global_n);
    ARRAY<T> local_sum(num_domains);
    
    // build interior views of x,b,p,z,temp
    threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,VECTOR_ND<T>&,VECTOR_ND<T>&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_One,domain_index,all_interior_indices,x,b,z_interior,x_interior,b_interior,p_interior,temp_interior);
    // adjust x for the null space
    if(enforce_compatibility && remove_null_space_solution_component){T sum=0;
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Sum,domain_index,all_interior_indices,x_interior,local_sum);
        for(int i=1;i<=num_domains;i++) sum+=local_sum(i);
        //threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,T>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Distribute,domain_index,all_interior_indices,x_interior,sum/global_n);}
        for(int i=1;i<=num_domains;i++) x_interior(i)-=sum/global_n;}
    
    // find initial residual, r=b-Ax - reusing b for the residual
    threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,const ARRAY<ARRAY<INTERVAL<int> > >&,SPARSE_MATRIX_FLAT_NXN<T>&,VECTOR_ND<T>&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_Two,domain_index,all_interior_indices,all_ghost_indices,A,x,b_interior,temp_interior);
    if(enforce_compatibility){T sum=0;
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Sum,domain_index,all_interior_indices,b_interior,local_sum);
        for(int i=1;i<=num_domains;i++) sum+=local_sum(i);
        //threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,T>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Distribute,domain_index,all_interior_indices,b_interior,sum/global_n);}
        for(int i=1;i<=num_domains;i++) b_interior(i)-=sum/global_n;}

    threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Max,domain_index,all_interior_indices,b_interior,local_sum);
    T max_val=0;for(int i=1;i<=num_domains;i++) max_val=max(max_val,local_sum(i));
    if(max_val<=global_tolerance){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(show_results) LOG::cout<<"NO ITERATIONS NEEDED"<<std::endl;
#endif
        return;}

    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
    ARRAY<SPARSE_MATRIX_FLAT_NXN<T>*> C(num_domains);
    threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,SPARSE_MATRIX_FLAT_NXN<T>&,ARRAY<SPARSE_MATRIX_FLAT_NXN<T>*>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_Three,domain_index,all_interior_indices,A,C);

    double rho=0,rho_old=0;
    for(int iteration=1;;iteration++){
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<SPARSE_MATRIX_FLAT_NXN<T>*>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_Four,domain_index,all_interior_indices,z_interior,b_interior,temp_interior,C);
        
        // for Neumann boundary conditions only, make sure z sums to zero
        if(enforce_compatibility){T sum=0;
            threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Sum,domain_index,all_interior_indices,z_interior,local_sum);
            for(int i=1;i<=num_domains;i++) sum+=local_sum(i);
            //threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,T>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Distribute,domain_index,all_interior_indices,z_interior,sum/global_n);}
            for(int i=1;i<=num_domains;i++) z_interior(i)-=sum/global_n;}
        
        // update search direction
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Dot,domain_index,all_interior_indices,z_interior,b_interior,local_sum);
        rho_old=rho;rho=0;for(int i=1;i<=num_domains;i++) rho+=local_sum(i);
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&,T,T,int>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_Five,domain_index,all_interior_indices,z_interior,p_interior,rho,rho_old,iteration);
        
        // update solution and residual
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,const ARRAY<ARRAY<INTERVAL<int> > >&,const SPARSE_MATRIX_FLAT_NXN<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_Six,domain_index,all_interior_indices,all_ghost_indices,A);
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Dot,domain_index,all_interior_indices,p_interior,temp_interior,local_sum);
        T sum=0;for(int i=1;i<=num_domains;i++) sum+=local_sum(i);
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<VECTOR_ND<T> >&,T>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Part_Seven,domain_index,all_interior_indices,x_interior,b_interior,p_interior,temp_interior,rho/sum);

        // remove null space component of b before computing residual norm because we might have converged up to the null space but have some null space component left due to roundoff
        if(enforce_compatibility){T sum=0;
            threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Sum,domain_index,all_interior_indices,b_interior,local_sum);
            for(int i=1;i<=num_domains;i++) sum+=local_sum(i);
            //threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,T>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Distribute,domain_index,all_interior_indices,b_interior,sum/global_n);}
            for(int i=1;i<=num_domains;i++) b_interior(i)-=sum/global_n;}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,ARRAY<VECTOR_ND<T> >&,ARRAY<T>&>(*this,&PCG_SPARSE_THREADED<TV>::Solve_Max,domain_index,all_interior_indices,b_interior,local_sum);
        T residual=0;for(int i=1;i<=num_domains;i++) residual=max(residual,local_sum(i));

        // check for convergence
        if(show_residual) LOG::cout<<residual<<std::endl;
        if(residual<=global_tolerance){if(show_results) LOG::cout<<"NUMBER OF ITERATIONS = "<<iteration<<std::endl;break;}
        if(iteration==desired_iterations){if(show_results) LOG::cout<<"DID NOT CONVERGE IN "<<iteration<<" ITERATIONS"<<std::endl;break;}
#endif
    }
    for(int i=1;i<=num_domains;i++) delete C(i);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_One(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,VECTOR_ND<T>& x,VECTOR_ND<T>& b,ARRAY<VECTOR_ND<T> >& z_interior,ARRAY<VECTOR_ND<T> >& x_interior,ARRAY<VECTOR_ND<T> >& b_interior,ARRAY<VECTOR_ND<T> >& p_interior,ARRAY<VECTOR_ND<T> >& temp_interior)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    assert(tid==domain_index(interior_domain.max_corner));
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);
    int interior_n=interior_indices.Size()+1;

    z_interior(tid).Resize(interior_n);
    x_interior(tid).Set_Subvector_View(x,interior_indices);
    b_interior(tid).Set_Subvector_View(b,interior_indices);
    p_interior(tid).Set_Subvector_View(p,interior_indices);
    temp_interior(tid).Set_Subvector_View(temp,interior_indices);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_Two(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,ARRAY<VECTOR_ND<T> >& b_interior,ARRAY<VECTOR_ND<T> >& temp_interior)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);
    const ARRAY<INTERVAL<int> >& ghost_indices=all_ghost_indices(tid);
 
    A.Times(interior_indices,ghost_indices,x,temp);
    b_interior(tid)-=temp_interior(tid);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_Three(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,ARRAY<SPARSE_MATRIX_FLAT_NXN<T>*>& C)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);
    
    C(tid)=0;
    if(incomplete_cholesky){
        C(tid)=A.Create_Submatrix(interior_indices);
        C(tid)->In_Place_Incomplete_Cholesky_Factorization(modified_incomplete_cholesky,modified_incomplete_cholesky_coefficient,
            preconditioner_zero_tolerance,preconditioner_zero_replacement);}
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_Four(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& z_interior,ARRAY<VECTOR_ND<T> >& b_interior,ARRAY<VECTOR_ND<T> >& temp_interior,ARRAY<SPARSE_MATRIX_FLAT_NXN<T>*>& C)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    
    if(incomplete_cholesky){
        // solve Mz=r
        C(tid)->Solve_Forward_Substitution(b_interior(tid),temp_interior(tid),true); // diagonal should be treated as the identity
        C(tid)->Solve_Backward_Substitution(temp_interior(tid),z_interior(tid),false,true);} // diagonal is inverted to save on divides
    else z_interior(tid)=b_interior(tid); // set z=r when there is no preconditioner
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_Five(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& z_interior,ARRAY<VECTOR_ND<T> >& p_interior,T rho,T rho_old,int iteration)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);
    int interior_n=interior_indices.Size()+1;
    
    T beta=0;if(iteration==1) p_interior(tid)=z_interior(tid);else{beta=(T)(rho/rho_old);for(int i=1;i<=interior_n;i++) p_interior(tid)(i)=z_interior(tid)(i)+beta*p_interior(tid)(i);}
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_Six(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,const SPARSE_MATRIX_FLAT_NXN<T>& A)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);
    const ARRAY<INTERVAL<int> >& ghost_indices=all_ghost_indices(tid);

    A.Times(interior_indices,ghost_indices,p,temp);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Part_Seven(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& x_interior,ARRAY<VECTOR_ND<T> >& b_interior,ARRAY<VECTOR_ND<T> >& p_interior,ARRAY<VECTOR_ND<T> >& temp_interior,T alpha)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);
    int interior_n=interior_indices.Size()+1;
   
    for(int i=1;i<=interior_n;i++){x_interior(tid)(i)+=alpha*p_interior(tid)(i);b_interior(tid)(i)-=alpha*temp_interior(tid)(i);}
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Distribute(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& interior,const T sum)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    interior(tid)-=sum;
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Sum(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& interior,ARRAY<T>& sum)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    sum(tid)=interior(tid).Sum_Double_Precision();
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Max(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& interior,ARRAY<T>& sum)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    sum(tid)=interior(tid).Max_Abs();
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Solve_Dot(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,ARRAY<VECTOR_ND<T> >& interior_1,ARRAY<VECTOR_ND<T> >& interior_2,ARRAY<T>& sum)
{
    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    sum(tid)=VECTOR_ND<T>::Dot_Product_Double_Precision(interior_1(tid),interior_2(tid));
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Subtract(VECTOR_ND<T>& vector,const T sum,int start_index,int end_index)
{
    for(int i=start_index;i<=end_index;i++) vector(i)-=sum;
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Sum(VECTOR_ND<T>& vector,ARRAY<T>& sum,int start_index,int end_index,int tid)
{
    sum(tid)=vector.Sum_Double_Precision(start_index,end_index);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Dot(VECTOR_ND<T>& vector1,VECTOR_ND<T>& vector2,ARRAY<T>& sum,int start_index,int end_index,int tid)
{
    sum(tid)=VECTOR_ND<T>::Dot_Product_Double_Precision(vector1,vector2,start_index,end_index);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Max(VECTOR_ND<T>& vector,ARRAY<T>& sum,int start_index,int end_index,int tid)
{
    sum(tid)=vector.Max_Abs(start_index,end_index);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Part_One(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T sum,int start_index,int end_index)
{
    for(int i=start_index;i<=end_index;i++) x(i)-=sum;
    A.Times(start_index,end_index,p,temp);
    for(int i=start_index;i<=end_index;i++) b(i)-=temp(i);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Part_Two(VECTOR_ND<T>& z,T rho,T rho_old,int iteration,int start_index,int end_index)
{
    T beta=0;if(iteration==1){for(int i=start_index;i<=end_index;i++) p(i)=z(i);}else{beta=(T)(rho/rho_old);for(int i=start_index;i<=end_index;i++) p(i)=z(i)+beta*p(i);}
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Part_Three(SPARSE_MATRIX_FLAT_NXN<T>& A,int start_index,int end_index)
{
    A.Times(start_index,end_index,p,temp);
}
template<class TV> void PCG_SPARSE_THREADED<TV>::
Threaded_Part_Four(VECTOR_ND<T>& x,VECTOR_ND<T>& b,T alpha,int start_index,int end_index)
{
    for(int i=start_index;i<=end_index;i++){x(i)+=alpha*p(i);b(i)-=alpha*temp(i);}
}
//#####################################################################
template class PCG_SPARSE_THREADED<VECTOR<float,1> >;
template class PCG_SPARSE_THREADED<VECTOR<float,2> >;
template class PCG_SPARSE_THREADED<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PCG_SPARSE_THREADED<VECTOR<double,1> >;
template class PCG_SPARSE_THREADED<VECTOR<double,2> >;
template class PCG_SPARSE_THREADED<VECTOR<double,3> >;
#endif
