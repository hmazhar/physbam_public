//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_ARRAY_ADAPTIVE
//#####################################################################
#ifndef __FACE_ARRAY_ADAPTIVE__
#define __FACE_ARRAY_ADAPTIVE__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
namespace PhysBAM{

template<class T,int d>
class FACE_ARRAY_ADAPTIVE:public ARRAY<T,FACE_INDEX<d> >
{
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;
    typedef ARRAY<T,FACE_INDEX<d> > BASE;
public:
    using BASE::Resize;

    ARRAY<ARRAY<T,FACE_INDEX<d> >*,TV_INT> sub_arrays;
    
    FACE_ARRAY_ADAPTIVE()
        :BASE()
    {}

    template<class T2>
    FACE_ARRAY_ADAPTIVE(const GRID<VECTOR<T2,d> >& grid,const int ghost_cells=0,const bool initialize_using_default_constructor=true)
        :BASE(grid,ghost_cells,initialize_using_default_constructor)
    {}

    FACE_ARRAY_ADAPTIVE(const RANGE<TV_INT>& domain_indices_input,const bool initialize_using_default_constructor=true)
        :BASE(domain_indices_input,initialize_using_default_constructor)
    {}

    template<class T2>
    FACE_ARRAY_ADAPTIVE(const ARRAY<T2,FACE_INDEX<d> >& old_array)
        :BASE(old_array)
    {}

    void Resize(const RANGE<TV_INT>& domain,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {BASE::Resize(domain,initialize_new_elements,copy_existing_elements,initialization_value);Initialize_Sub_Arrays(domain);}

    void Initialize_Sub_Arrays(const RANGE<TV_INT>& domain)
    {sub_arrays.Resize(domain);}

//#####################################################################
};
}
#endif
