// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#pragma once

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

using namespace PhysBAM;

namespace PhysBAM{

template<class T>
class CG_VECTOR:public KRYLOV_VECTOR_BASE<T>
{
    typedef KRYLOV_VECTOR_BASE<T> BASE;

    ARRAY<VECTOR<T,3> >& array;
public:
    CG_VECTOR(ARRAY<VECTOR<T,3> >& array_input) : array(array_input) {}

    static ARRAY<VECTOR<T,3> >& Array(BASE& base_array)
    {return ((CG_VECTOR&)(base_array)).array;}

    static const ARRAY<VECTOR<T,3> >& Array(const BASE& base_array)
    {return ((const CG_VECTOR&)(base_array)).array;}

    BASE& operator+=(const BASE& bv)
    {array+=Array(bv);return *this;}

    BASE& operator-=(const BASE& bv)
    {array-=Array(bv);return *this;}

    BASE& operator*=(const T a)
    {array*=a;return *this;}

    void Copy(const T c,const BASE& bv)
    {ARRAY<VECTOR<T,3> >::Copy(c,Array(bv),array);}

    void Copy(const T c1,const BASE& bv1,const BASE& bv2)
    {ARRAY<VECTOR<T,3> >::Copy(c1,Array(bv1),Array(bv2),array);}

    int Raw_Size() const
    {return array.Flattened().m;}
    
    T& Raw_Get(int i)
    {return array.Flattened()(i);}
};

}
