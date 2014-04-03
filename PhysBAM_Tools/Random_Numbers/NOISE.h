//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NOISE
//##################################################################### 
#ifndef __NOISE_h__
#define __NOISE_h__
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class NOISE 
{
    static const int noise_seed,prime1,prime2,prime3;
    static const T phi;
    static T points[];
    static bool initialized;

    NOISE()
    {}

//#####################################################################
    static int Hash(const int a,const int b,const int c);
    static void Initialize_Noise();
public:
    static T Noise_Function1(const VECTOR<T,3>& p);
    static VECTOR<T,3> Noise_Function3(const VECTOR<T,3>& p);
    static T Noise1(const VECTOR<T,3>& p,const int octaves,const T persistance=.5); // TODO: Remove these old versions
    static void Noise3(const VECTOR<T,3>& p,VECTOR<T,3>& v,const int octaves,const T persistance=.5); // TODO: Remove these old versions
    static T Normalization_Factor(const int octaves,const T persistance);
    static T Fractional_Brownian(const VECTOR<T,3>& p,const int octaves,const T lacunarity,const T gain);
    static T Turbulence(const VECTOR<T,3>& p,const int octaves,const T lacunarity,const T gain);
    static VECTOR<T,3> Fractional_Brownian_3D(const VECTOR<T,3>& p,const int octaves,const T lacunarity,const T gain);
    static VECTOR<T,3> Turbulence_3D(const VECTOR<T,3>& p,const int octaves,const T lacunarity,const T gain);
//#####################################################################
};
}
#endif
