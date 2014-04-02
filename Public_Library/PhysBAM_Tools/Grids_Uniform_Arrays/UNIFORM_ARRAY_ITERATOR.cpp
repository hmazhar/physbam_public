//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_ARRAY_ITERATOR
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/UNIFORM_ARRAY_ITERATOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<int dimension> UNIFORM_ARRAY_ITERATOR<dimension>::
UNIFORM_ARRAY_ITERATOR(const RANGE<TV_INT>& domain_input)
    :domain(domain_input)
{
    Reset();
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<int dimension> void UNIFORM_ARRAY_ITERATOR<dimension>::
Next_Helper()
{
    index(dimension)=domain.min_corner(dimension);
    for(int i=dimension-1;i>=1;i--){
        if(index(i)<domain.max_corner(i)){index(i)++;return;}
        index(i)=domain.min_corner(i);}
    valid=false;
}
//#####################################################################
template class UNIFORM_ARRAY_ITERATOR<0>;
template class UNIFORM_ARRAY_ITERATOR<1>;
template class UNIFORM_ARRAY_ITERATOR<2>;
template class UNIFORM_ARRAY_ITERATOR<3>;
