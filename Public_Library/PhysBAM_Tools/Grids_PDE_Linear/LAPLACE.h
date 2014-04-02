//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE
//#####################################################################
#ifndef __LAPLACE__
#define __LAPLACE__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
namespace PhysBAM{

template<class T>
class LAPLACE:public NONCOPYABLE
{
public:
    T tolerance;
    bool solve_neumann_regions;
    int number_of_regions;
    T relative_tolerance;
    T absolute_tolerance;

    LAPLACE()
        :number_of_regions(0)
    {
        Set_Relative_Tolerance();
        Set_Absolute_Tolerance();
        Solve_Neumann_Regions();
    }

    virtual ~LAPLACE()
    {}

    void Set_Relative_Tolerance(const T relative_tolerance_input=1e-7)
    {relative_tolerance=relative_tolerance_input;}

    void Set_Absolute_Tolerance(const T absolute_tolerance_input=1e-14)
    {absolute_tolerance=absolute_tolerance_input;}

    void Solve_Neumann_Regions(const bool solve_neumann_regions_input=true)
    {solve_neumann_regions=solve_neumann_regions_input;}

    virtual void Compute_beta_And_Add_Jumps_To_b(const T dt,const T time)
    {}

    void Find_Tolerance(const VECTOR_ND<T>& b)
    {tolerance=max(relative_tolerance*b.Max_Abs(),absolute_tolerance);}

//#####################################################################
};
}
#endif
