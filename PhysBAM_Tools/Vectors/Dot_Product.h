//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Dot_Product 
//#####################################################################
#ifndef __Dot_Product__
#define __Dot_Product__

namespace PhysBAM{

template<class TV>
inline typename TV::SCALAR Dot_Product(const TV& v1,const TV& v2)
{return TV::Dot_Product(v1,v2);}

inline float Dot_Product(const float a1,const float a2)
{return a1*a2;}

inline double Dot_Product(const double a1,const double a2)
{return a1*a2;}

template<class TV>
inline double Dot_Product_Double_Precision(const TV& v1,const TV& v2)
{return TV::Dot_Product_Double_Precision(v1,v2);}

template<class T,int d> class VECTOR;
template<class T,int d>
inline double Dot_Product_Double_Precision(const VECTOR<T,d>& v1,const VECTOR<T,d>& v2)
{return VECTOR<T,d>::Dot_Product(v1,v2);}

inline double Dot_Product_Double_Precision(const float a1,const float a2)
{return a1*a2;}

inline double Dot_Product_Double_Precision(const double a1,const double a2)
{return a1*a2;}

}
#endif
