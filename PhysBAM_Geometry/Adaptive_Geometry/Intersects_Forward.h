//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INTERSECTS_FORWARD__
#define __INTERSECTS_FORWARD__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{
namespace GEOMETRIC_PREDICATES_DETAIL{
//#####################################################################
// Intersects
//#####################################################################
template<class EXACT_TYPE,class T,int D>
bool Intersects(const VECTOR<T,D>& point,const VECTOR<VECTOR<T,D>,D+1>& simplex,bool* is_degenerate_p=0);

template<class EXACT_TYPE,class T,int D>
bool Intersects(const VECTOR<VECTOR<T,D>,D+1>& simplex,const VECTOR<T,D>& point,bool* is_degenerate_p=0);

template<class EXACT_TYPE,class T,int D,int N1,int N2>
typename ENABLE_IF<(N1<N2),bool>::TYPE
Intersects(const VECTOR<VECTOR<T,D>,N1>& simplex1,const VECTOR<VECTOR<T,D>,N2>& simplex2,bool* is_degenerate_p=0);

template<class EXACT_TYPE,class T,int D>
bool Intersects(const VECTOR<VECTOR<T,D>,D+1>& simplex1,const VECTOR<VECTOR<T,D>,1>& point,bool* is_degenerate_p=0);

template<class EXACT_TYPE,class T,int D,int N1,int N2>
typename ENABLE_IF<(N2<=N1&&N1<=D&&N1+N2<=D+2),bool>::TYPE
Intersects(const VECTOR<VECTOR<T,D>,N1>& simplex1,const VECTOR<VECTOR<T,D>,N2>& simplex2,bool* is_degenerate_p);

template<class EXACT_TYPE,class T,int D,int N1,int N2>
typename ENABLE_IF<(N2<=N1&&N1<=D&&N1+N2>=D+3),bool>::TYPE
Intersects(const VECTOR<VECTOR<T,D>,N1>& simplex1,const VECTOR<VECTOR<T,D>,N2>& simplex2,bool* is_degenerate_p);

template<class EXACT_TYPE,class T,int D,int N2>
bool Intersects(const VECTOR<VECTOR<T,D>,D+1>& simplex1,const VECTOR<VECTOR<T,D>,N2>& simplex2,bool* is_degenerate_p=0);

template<class EXACT_TYPE,class T>
bool Intersects(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,bool* is_degenerate_p=0);
//#####################################################################
}
}
#endif
