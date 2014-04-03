//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_ARRAY_ITERATOR
//#####################################################################
#ifndef __UNIFORM_ARRAY_ITERATOR__
#define __UNIFORM_ARRAY_ITERATOR__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<int dimension>
class UNIFORM_ARRAY_ITERATOR
{
    typedef VECTOR<int,dimension> TV_INT;
public:
    const RANGE<TV_INT>& domain;
    TV_INT index;
    bool valid;

    UNIFORM_ARRAY_ITERATOR(const RANGE<TV_INT>& domain_input);

    void Reset()
    {if(dimension==0){valid=true;index=TV_INT();return;}
    valid=true;index=domain.Minimum_Corner();}

    bool Valid() const
    {return valid;}

    void Next() PHYSBAM_ALWAYS_INLINE
    {if(dimension==0){valid=false;return;}
    if(index(dimension)<domain.max_corner(dimension)) index(dimension)++;else Next_Helper();}
    
    const TV_INT& Index()
    {return index;}

protected:
    void Next_Helper();
};
//#####################################################################
}
#endif
