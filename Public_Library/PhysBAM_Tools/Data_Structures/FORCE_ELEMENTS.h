//#####################################################################
// Copyright 2007-2009, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FORCE_ELEMENTS
//#####################################################################
#ifndef __FORCE_ELEMENTS__
#define __FORCE_ELEMENTS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>

namespace PhysBAM{
template<class TV> class MPI_SOLIDS;
class FORCE_ELEMENTS
{
public:
    ARRAY<int> indices;

    FORCE_ELEMENTS();
    ~FORCE_ELEMENTS();

    template<int d,class T_ARRAY> void Update(const ARRAY_BASE<VECTOR<int,d>,T_ARRAY>& base_elements,const ARRAY<bool>& particle_is_simulated);
    template<int d,class T_ARRAY> void Update(const ARRAY_BASE<VECTOR<int,d>,T_ARRAY>& base_elements);
    template<class T_ARRAY> void Update(const ARRAY_BASE<int,T_ARRAY>& base_elements,const ARRAY<bool>& particle_is_simulated);
    template<class T_ARRAY> void Update(const ARRAY_BASE<int,T_ARRAY>& base_elements);
    void Print() const;

    class ITERATOR
    {
    public:
        const ARRAY<int>& indices;
        int index;

        ITERATOR(const FORCE_ELEMENTS& force_elements)
            :indices(force_elements.indices),index(1)
        {}

        void Next()
        {index++;}

        bool Valid() const
        {return index<=indices.m;}

        int Data() const
        {return indices(index);}
    };
};
}
#endif
