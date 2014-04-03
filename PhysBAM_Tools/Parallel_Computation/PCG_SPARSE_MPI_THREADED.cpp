//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifdef USE_MPI
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/INT_ITERATOR_THREADED.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI_THREADED.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>
using namespace PhysBAM;
//#####################################################################
// Function Parallel_Solve
//#####################################################################
template<class T> void PCG_SPARSE_MPI_THREADED<T>::
Solve(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,const ARRAY<INTERVAL<int> >& all_interior_indices,const ARRAY<ARRAY<INTERVAL<int> > >& all_ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance)
{
    pcg_threaded.Init_Barriers();
    BASE::Initialize_Datatypes();

    RANGE<TV_INT> interior_domain(domain);interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
    int tid=domain_index(interior_domain.min_corner);
    assert(tid==domain_index(interior_domain.max_corner));
    const INTERVAL<int>& interior_indices=all_interior_indices(tid);
    const ARRAY<INTERVAL<int> >& ghost_indices=all_ghost_indices(tid);

    int local_n=A.n,interior_n=interior_indices.Size()+1;
    int global_n=Global_Sum_Int(interior_n);
    T global_tolerance=Global_Max(tolerance);
    pcg.maximum_iterations=40;
    int desired_iterations=global_n;if(pcg.enforce_compatibility) desired_iterations--;if(pcg.maximum_iterations) desired_iterations=min(desired_iterations,pcg.maximum_iterations);

    VECTOR_ND<T>& temp=pcg_threaded.temp;
    VECTOR_ND<T>& p=pcg_threaded.p;
    VECTOR_ND<T> z_interior(interior_n,false);
#ifdef USE_PTHREADS
    pthread_mutex_lock(&pcg_threaded.lock);
#endif
    temp.Resize(local_n);p.Resize(local_n);
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&pcg_threaded.lock);
#endif

    // build interior views of x,b,p,z,temp
    VECTOR_ND<T> x_interior,b_interior,p_interior,temp_interior;
    x_interior.Set_Subvector_View(x,interior_indices);
    b_interior.Set_Subvector_View(b,interior_indices);
    p_interior.Set_Subvector_View(p,interior_indices);
    temp_interior.Set_Subvector_View(temp,interior_indices);

    // adjust x for the null space
    if(pcg.enforce_compatibility && pcg.remove_null_space_solution_component) x_interior-=(T)(Global_Sum(x_interior.Sum_Double_Precision(),tid)/global_n);

    // find initial residual, r=b-Ax - reusing b for the residual
    Fill_Ghost_Cells(x);
    A.Times(interior_indices,ghost_indices,x,temp);
    b_interior-=temp_interior;
    if(pcg.enforce_compatibility) b_interior-=(T)(Global_Sum(b_interior.Sum_Double_Precision(),tid)/global_n);
    if(Global_Max(b_interior.Max_Abs())<=global_tolerance){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(pcg.show_results) LOG::cout<<"NO ITERATIONS NEEDED"<<std::endl;
#endif
        return;}

    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
    SPARSE_MATRIX_FLAT_NXN<T>* C=0;
    if(pcg.incomplete_cholesky){
        C=A.Create_Submatrix(interior_indices);
        C->In_Place_Incomplete_Cholesky_Factorization(pcg.modified_incomplete_cholesky,pcg.modified_incomplete_cholesky_coefficient,
            pcg.preconditioner_zero_tolerance,pcg.preconditioner_zero_replacement);}

    double rho=0,rho_old=0;
    for(int iteration=1;;iteration++){
        //if(show_results) LOG::Time("Iteration");
        if(pcg.incomplete_cholesky){
            // solve Mz=r
            C->Solve_Forward_Substitution(b_interior,temp_interior,true); // diagonal should be treated as the identity
            C->Solve_Backward_Substitution(temp_interior,z_interior,false,true);} // diagonal is inverted to save on divides
        else z_interior=b_interior; // set z=r when there is no preconditioner

        // for Neumann boundary conditions only, make sure z sums to zero
        if(pcg.enforce_compatibility) z_interior-=(T)(Global_Sum(z_interior.Sum_Double_Precision(),tid)/global_n);

        // update search direction
        rho_old=rho;rho=Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(z_interior,b_interior),tid);
        T beta=0;if(iteration==1) p_interior=z_interior;else{beta=(T)(rho/rho_old);for(int i=1;i<=interior_n;i++) p_interior(i)=z_interior(i)+beta*p_interior(i);} // when iteration=1, beta=0

        // update solution and residual
        Fill_Ghost_Cells(p);
        A.Times(interior_indices,ghost_indices,p,temp);
        T alpha=(T)(rho/Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(p_interior,temp_interior),tid));
        for(int i=1;i<=interior_n;i++){x_interior(i)+=alpha*p_interior(i);b_interior(i)-=alpha*temp_interior(i);}

        // remove null space component of b before computing residual norm because we might have converged up to the null space but have some null space component left due to roundoff
        if(pcg.enforce_compatibility) b_interior-=(T)(Global_Sum(b_interior.Sum_Double_Precision(),tid)/global_n);

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        T residual=Global_Max(b_interior.Max_Abs());

        // check for convergence
        if(pcg.show_residual) LOG::cout<<residual<<std::endl;
        if(residual<=global_tolerance){if(pcg.show_results) LOG::cout<<"NUMBER OF ITERATIONS = "<<iteration<<std::endl;break;}
        if(iteration==desired_iterations){if(pcg.show_results) LOG::cout<<"DID NOT CONVERGE IN "<<iteration<<" ITERATIONS"<<std::endl;break;}
#endif
    }
    if(pcg.show_results) LOG::Time("Done");
    Fill_Ghost_Cells(x);
}
//#####################################################################
template class PCG_SPARSE_MPI_THREADED<VECTOR<float,1> >;
template class PCG_SPARSE_MPI_THREADED<VECTOR<float,2> >;
template class PCG_SPARSE_MPI_THREADED<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PCG_SPARSE_MPI_THREADED<VECTOR<double,1> >;
template class PCG_SPARSE_MPI_THREADED<VECTOR<double,2> >;
template class PCG_SPARSE_MPI_THREADED<VECTOR<double,3> >;
#endif
#endif
