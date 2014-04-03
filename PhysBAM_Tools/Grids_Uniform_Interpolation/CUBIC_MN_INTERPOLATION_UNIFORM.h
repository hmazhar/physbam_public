//#####################################################################
// Copyright 2003-2010, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUBIC_MN_INTERPOLATION_UNIFORM 
//#####################################################################
#ifndef __CUBIC_MN_INTERPOLATION_UNIFORM__
#define __CUBIC_MN_INTERPOLATION_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Interpolation/CUBIC_MN_INTERPOLATION.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<T_GRID>
class CUBIC_MN_INTERPOLATION_UNIFORM:public INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    //typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> > T_ARRAYS_T2;
public:
    typedef INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP> BASE;
    using BASE::Clamped_Index_Interior_End_Minus_One;

    T b,c; // two parameter family
    CUBIC_MN_INTERPOLATION<T,T2> cubic_mn_interpolation;

    CUBIC_MN_INTERPOLATION_UNIFORM();
    ~CUBIC_MN_INTERPOLATION_UNIFORM();

    void Set_Sharpness(T sharpness_input=2./3)
    {Set_Parameters(1-sharpness_input,(T).5*sharpness_input);}

    void Set_Parameters(const T b_input=1./3,const T c_input=1./3)
    {b=b_input;c=c_input;cubic_mn_interpolation.Set_Parameters(b,c);}

    T2 Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE;
    ARRAY<PAIR<TV_INT,T> > Clamped_To_Array_Weights(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE;
    T2 Periodic(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const;
    // T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const;
    // T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const;
    // T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const;
    // ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const;
    // ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const;
    // ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const;
    // only works for power of two grids!
//    T2 From_Base_Node_Periodic(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<T2,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const;
//#####################################################################
};
}
#endif
