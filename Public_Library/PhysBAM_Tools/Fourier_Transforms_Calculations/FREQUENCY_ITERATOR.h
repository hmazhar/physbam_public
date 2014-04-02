//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FREQUENCY_ITERATOR
//#####################################################################
#ifndef __FREQUENCY_ITERATOR__
#define __FREQUENCY_ITERATOR__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class FREQUENCY_ITERATOR:public GRID<TV>::BASE_ITERATOR
{
private:
    typedef typename TV::SCALAR T;typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename GRID<TV>::BASE_ITERATOR T_ITERATOR;

    using T_ITERATOR::grid;using T_ITERATOR::index;
    const TV coefficients;
public:
    struct HALF_DOMAIN{};struct FULL_DOMAIN{};

    FREQUENCY_ITERATOR(const GRID<TV>& grid,HALF_DOMAIN domain_type=HALF_DOMAIN())
        :T_ITERATOR(grid,Domain_Indices(grid,domain_type)),coefficients((T)(2*pi)/grid.domain.Edge_Lengths())
    {}

    FREQUENCY_ITERATOR(const GRID<TV>& grid,FULL_DOMAIN domain_type)
        :T_ITERATOR(grid,Domain_Indices(grid,domain_type)),coefficients((T)(2*pi)/grid.domain.Edge_Lengths())
    {}

    static RANGE<VECTOR<int,1> > Domain_Indices(const GRID<VECTOR<T,1> >& grid,HALF_DOMAIN domain_type=HALF_DOMAIN())
    {return RANGE<VECTOR<int,1> >(0,grid.counts.x/2);}

    static RANGE<VECTOR<int,1> > Domain_Indices(const GRID<VECTOR<T,1> >& grid,FULL_DOMAIN)
    {return RANGE<VECTOR<int,1> >(0,grid.counts.x);}

    static RANGE<VECTOR<int,2> > Domain_Indices(const GRID<VECTOR<T,2> >& grid,HALF_DOMAIN domain_type=HALF_DOMAIN())
    {return RANGE<VECTOR<int,2> >(0,grid.counts.x-1,0,grid.counts.y/2);}

    static RANGE<VECTOR<int,2> > Domain_Indices(const GRID<VECTOR<T,2> >& grid,FULL_DOMAIN)
    {return RANGE<VECTOR<int,2> >(0,grid.counts.x,0,grid.counts.y);}

    TV_INT Index() const
    {return index;}

    TV Frequency() const
    {TV_INT counts=grid.counts;TV k;
    for(int i=1;i<=k.m;i++) k[i]=coefficients[i]*(index[i]<=counts[i]/2?index[i]:index[i]-counts[i]);
    return k;}
};
//#####################################################################
}
#endif
