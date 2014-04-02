//#####################################################################
// Copyright 2004-2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRAL_GROUP
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __TETRAHEDRAL_GROUP__
#define __TETRAHEDRAL_GROUP__

#include <PhysBAM_Tools/Matrices/ROTATION.h>
namespace PhysBAM{

template<class T>
class TETRAHEDRAL_GROUP
{
public:
    int g; // group element (0 to 11)
private:
    static const ROTATION<VECTOR<T,3> > quaternions[12];
    static const ROTATION<VECTOR<int,3> > quaternions_times_2[12];
    static const int multiplication_table[12][12];
    static const int inversion_table[12];
public:
    static const TETRAHEDRAL_GROUP<T> e,i,j,k,a,a2,b,b2,c,c2,d,d2; // letter names for elements

    TETRAHEDRAL_GROUP(const int g_input=0)
        :g(g_input)
    {
        assert(0<=g && g<12);
    }

    TETRAHEDRAL_GROUP<T> operator*(const TETRAHEDRAL_GROUP<T> g2) const
    {return multiplication_table[g][g2.g];}

    TETRAHEDRAL_GROUP<T> Inverse() const
    {return inversion_table[g];}

    VECTOR<T,3> Rotate(const VECTOR<T,3>& v) const
    {return quaternions[g].Rotate(v);}

    VECTOR<int,3> Rotate(const VECTOR<int,3>& v) const
    {return quaternions_times_2[g].Rotate(v)/4;}

    VECTOR<T,3> Inverse_Rotate(const VECTOR<T,3>& v) const
    {return quaternions[g].Inverse_Rotate(v);}

    VECTOR<int,3> Inverse_Rotate(const VECTOR<int,3>& v) const
    {return quaternions_times_2[g].Inverse_Rotate(v)/4;}

    static TETRAHEDRAL_GROUP<T> Cyclic_Shift_Axes(const int times=1)
    {switch(times%3){case 0:return e;case 1:return a2;default:return a;}}

    static void Assert_Correctness()
    {for(int g=0;g<12;g++)assert(!multiplication_table[g][inversion_table[g]] && !multiplication_table[inversion_table[g]][g]);
    for(int g=0;g<12;g++)for(int h=0;h<12;h++){
        ROTATION<VECTOR<T,3> > gh=quaternions[multiplication_table[g][h]],g_times_h=quaternions[g]*quaternions[h];
        assert(gh==g_times_h);}}

//#####################################################################
};
}
#endif

#endif
