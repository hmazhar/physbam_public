//#####################################################################
// Copyright 2007-2009, Michael Lentine, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
FORCE_ELEMENTS::
FORCE_ELEMENTS()
{}
//#####################################################################
// Destructor
//#####################################################################
FORCE_ELEMENTS::
~FORCE_ELEMENTS()
{}
//#####################################################################
// Function Update
//#####################################################################
template<int d,class T_ARRAY> void FORCE_ELEMENTS::
Update(const ARRAY_BASE<VECTOR<int,d>,T_ARRAY>& base_elements)
{
    indices=IDENTITY_ARRAY<>(base_elements.Size());
}
//#####################################################################
// Function Update
//#####################################################################
template<int d,class T_ARRAY> void FORCE_ELEMENTS::
Update(const ARRAY_BASE<VECTOR<int,d>,T_ARRAY>& base_elements,const ARRAY<bool>& particle_is_simulated)
{
    const T_ARRAY& elements=base_elements.Derived();
    indices.Remove_All();
    for(int i=1;i<=elements.Size();i++)
        if(particle_is_simulated.Subset(elements(i)).Contains(true))
            indices.Append(i);
}
//#####################################################################
// Function Update
//#####################################################################
template<class T_ARRAY> void FORCE_ELEMENTS::
Update(const ARRAY_BASE<int,T_ARRAY>& base_elements,const ARRAY<bool>& particle_is_simulated)
{
    HASHTABLE<int> found;
    const T_ARRAY& elements=base_elements.Derived();
    indices.Remove_All();
    for(int i=1;i<=elements.Size();i++)
        if(particle_is_simulated(elements(i)))
            if(found.Set(elements(i)))
                indices.Append(elements(i));
}
//#####################################################################
// Function Update
//#####################################################################
template<class T_ARRAY> void FORCE_ELEMENTS::
Update(const ARRAY_BASE<int,T_ARRAY>& base_elements)
{
    HASHTABLE<int> found;
    const T_ARRAY& elements=base_elements.Derived();
    indices.Remove_All();
    for(int i=1;i<=elements.Size();i++)
        if(found.Set(elements(i)))
            indices.Append(elements(i));
}
//#####################################################################
// Function Print
//#####################################################################
void FORCE_ELEMENTS::
Print() const
{
    LOG::cout<<"indices: "<<indices<<std::endl;
}
//#####################################################################
#define INSTANTIATION_HELPER_d(d) \
    template void FORCE_ELEMENTS::Update<d,ARRAY<VECTOR<int,d> > >(ARRAY_BASE<VECTOR<int,d>,ARRAY<VECTOR<int,d> > > const&); \
    template void FORCE_ELEMENTS::Update<d,ARRAY<VECTOR<int,d> > >(ARRAY_BASE<VECTOR<int,d>,ARRAY<VECTOR<int,d> > > const&,ARRAY<bool> const&);
INSTANTIATION_HELPER_d(2);INSTANTIATION_HELPER_d(3);INSTANTIATION_HELPER_d(4);INSTANTIATION_HELPER_d(8);
#define INSTANTIATION_HELPER(T_ARRAY) \
    template void FORCE_ELEMENTS::Update<T_ARRAY<int> >(ARRAY_BASE<int,T_ARRAY<int> > const&); \
    template void FORCE_ELEMENTS::Update<T_ARRAY<int> >(ARRAY_BASE<int,T_ARRAY<int> > const&,ARRAY<bool> const&);
INSTANTIATION_HELPER(ARRAY);INSTANTIATION_HELPER(ARRAY_VIEW);INSTANTIATION_HELPER(IDENTITY_ARRAY);
