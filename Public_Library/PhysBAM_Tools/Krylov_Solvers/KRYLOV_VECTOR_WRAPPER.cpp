#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_WRAPPER<T,TV>::
KRYLOV_VECTOR_WRAPPER()
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_WRAPPER<T,TV>::
KRYLOV_VECTOR_WRAPPER(TV vector)
    :v(vector)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,class TV> template<class VECTOR,class INDICES> KRYLOV_VECTOR_WRAPPER<T,TV>::
KRYLOV_VECTOR_WRAPPER(VECTOR& vector,INDICES& index)
    :v(vector,index)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_WRAPPER<T,TV>::
~KRYLOV_VECTOR_WRAPPER()
{
}
//#####################################################################
// Operator +=
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_BASE<T>& KRYLOV_VECTOR_WRAPPER<T,TV>::
operator+=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v+=dynamic_cast<const KRYLOV_VECTOR_WRAPPER&>(bv).v;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_BASE<T>& KRYLOV_VECTOR_WRAPPER<T,TV>::
operator-=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v-=dynamic_cast<const KRYLOV_VECTOR_WRAPPER&>(bv).v;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_BASE<T>& KRYLOV_VECTOR_WRAPPER<T,TV>::
operator*=(const T a)
{
    v*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class T,class TV> void KRYLOV_VECTOR_WRAPPER<T,TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv)
{
    v.Copy(c,dynamic_cast<const KRYLOV_VECTOR_WRAPPER&>(bv).v,v);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class T,class TV> void KRYLOV_VECTOR_WRAPPER<T,TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2)
{
    v.Copy(c,dynamic_cast<const KRYLOV_VECTOR_WRAPPER&>(bv1).v,dynamic_cast<const KRYLOV_VECTOR_WRAPPER&>(bv2).v,v);
}
namespace{
inline int Raw_Size_Helper(const float& p){return 1;}
inline int Raw_Size_Helper(const double& p){return 1;}
template<class T,class T_VECTOR> inline int Raw_Size_Helper(const VECTOR_BASE<T,T_VECTOR>& p){if(!p.Size()) return 0;return p.Size()*Raw_Size_Helper(p(1));}
template<class T,class T_ARRAY> inline int Raw_Size_Helper(const ARRAY_BASE<T,T_ARRAY>& p){if(!p.Size()) return 0;return p.Size()*Raw_Size_Helper(p(1));}
inline float& Raw_Get_Helper(int i,float*,float& p){assert(i==1);return p;}
inline double& Raw_Get_Helper(int i,double*,double& p){assert(i==1);return p;}
template<class S,class T,class T_VECTOR> inline S& Raw_Get_Helper(int i,S* a,VECTOR_BASE<T,T_VECTOR>& p){int s=Raw_Size_Helper(p(1));return Raw_Get_Helper((i-1)%s+1,a,p((i-1)/s+1));}
template<class S,class T,class T_ARRAY> inline S& Raw_Get_Helper(int i,S* a,ARRAY_BASE<T,T_ARRAY>& p){int s=Raw_Size_Helper(p(1));return Raw_Get_Helper((i-1)%s+1,a,p((i-1)/s+1));}
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class T,class TV> int KRYLOV_VECTOR_WRAPPER<T,TV>::
Raw_Size() const
{
    return Raw_Size_Helper(v);
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class T,class TV> T& KRYLOV_VECTOR_WRAPPER<T,TV>::
Raw_Get(int i)
{
    return Raw_Get_Helper(i,(T*)0,v);
}
template class KRYLOV_VECTOR_WRAPPER<float,ARRAY<VECTOR_ND<float> > >;
template class KRYLOV_VECTOR_WRAPPER<float,VECTOR_ND<float> >;
template KRYLOV_VECTOR_BASE<float>& KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<float,int>,ARRAY<int>&> >::operator*=(float);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<VECTOR<float,2>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<float,2>,int>&,ARRAY<int> const&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<VECTOR<float,2>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<VECTOR<float,3>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<float,3>,int>&,ARRAY<int> const&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<VECTOR<float,3>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<float,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<float,int>&,ARRAY<int> const&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<float,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<float,int>&,ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<float,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,2>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<float,2>,int>&,ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,2>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY_VIEW<VECTOR<float,2>,int>&,ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,2>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,3>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<float,3>,int>&,ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,3>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY_VIEW<VECTOR<float,3>,int>&,ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,3>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<float,VECTOR_ND<float>&>::KRYLOV_VECTOR_WRAPPER(VECTOR_ND<float>&);
template KRYLOV_VECTOR_WRAPPER<float,VECTOR_ND<float>&>::~KRYLOV_VECTOR_WRAPPER();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class KRYLOV_VECTOR_WRAPPER<double,ARRAY<VECTOR_ND<double> > >;
template class KRYLOV_VECTOR_WRAPPER<double,VECTOR_ND<double> >;
template KRYLOV_VECTOR_BASE<double>& KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<double,int>,ARRAY<int>&> >::operator*=(double);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<VECTOR<double,2>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<double,2>,int>&,ARRAY<int> const&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<VECTOR<double,2>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<VECTOR<double,3>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<double,3>,int>&,ARRAY<int> const&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<VECTOR<double,3>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<double,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<double,int>&,ARRAY<int> const&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<double,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<double,int>&,ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<double,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,2>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<double,2>,int>&,ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,2>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY_VIEW<VECTOR<double,2>,int>&,ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,2>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,3>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<double,3>,int>&,ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,3>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY_VIEW<VECTOR<double,3>,int>&,ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,3>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<double,VECTOR_ND<double>&>::KRYLOV_VECTOR_WRAPPER(VECTOR_ND<double>&);
template KRYLOV_VECTOR_WRAPPER<double,VECTOR_ND<double>&>::~KRYLOV_VECTOR_WRAPPER();
#endif
