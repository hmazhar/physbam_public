//#####################################################################
// Copyright 2008, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __POLYGON_ORIENTATION__
#define __POLYGON_ORIENTATION__
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_ATOM.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_OBJECT.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_FLOAT.h>

namespace PhysBAM{
//##################################################################### 
// Function Polygon_Orientation
//##################################################################### 
template<class T_ARRAY_TV2>
int Polygon_Orientation(const T_ARRAY_TV2& polygon_coordinates)
{
    typedef typename T_ARRAY_TV2::ELEMENT TV2;
    typedef typename TV2::ELEMENT T;
    typedef EXACT_FLOAT<T> T_EXACT;
    typedef ADAPTIVE_ATOM<T,T_EXACT> ATOM_TYPE;
    if(polygon_coordinates.Size()<3) return 0;
    ADAPTIVE_OBJECT<T_EXACT> signed_area_x2(ATOM_TYPE(0.0));
    for(int p=1;p<=polygon_coordinates.Size();++p){
        int q=p%polygon_coordinates.Size()+1;
        const TV2 &xp=polygon_coordinates(p),&xq=polygon_coordinates(q);
        signed_area_x2+=(ATOM_TYPE(xp[1])*ATOM_TYPE(xq[2])-ATOM_TYPE(xp[2])*ATOM_TYPE(xq[1]));}
    return signed_area_x2.Sign();
//##################################################################### 
}
}
#endif

