//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Geoffrey Irving, Mike Rodgers, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Hash
//#####################################################################
#ifndef __Hash__
#define __Hash__

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/SIDED_FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <cstring>
#include <string>
namespace PhysBAM{

template<class TV> class RANGE;

//#####################################################################
// Class HASH
//#####################################################################
class HASH
{
    static const int missing_element_hash=32138912;
public:

    int value;

    explicit HASH()
        :value(missing_element_hash)
    {}

    explicit HASH(const int key)
        :value(int_hash(key))
    {}

    template<class T1,class T2>
    HASH(const T1& key1,const T2& key2)
        :value(triple_int_hash(missing_element_hash,Value(key1),Value(key2)))
    {}

    template<class T1,class T2,class T3>
    HASH(const T1& key1,const T2& key2,const T3& key3)
        :value(triple_int_hash(Value(key1),Value(key2),Value(key3)))
    {}

    template<class T1,class T2,class T3,class T4>
    HASH(const T1& key1,const T2& key2,const T3& key3,const T4& key4)
        :value(triple_int_hash(triple_int_hash(missing_element_hash,Value(key1),Value(key2)),Value(key3),Value(key4)))
    {}

    template<class T1,class T2,class T3,class T4,class T5>
    HASH(const T1& key1,const T2& key2,const T3& key3,const T4& key4,const T5& key5)
        :value(triple_int_hash(triple_int_hash(Value(key1),Value(key2),Value(key3)),Value(key4),Value(key5)))
    {}

    template<class T1,class T2,class T3,class T4,class T5,class T6>
    HASH(const T1& key1,const T2& key2,const T3& key3,const T4& key4,const T5& key5,const T6& key6)
        :value(triple_int_hash(triple_int_hash(triple_int_hash(missing_element_hash,Value(key1),Value(key2)),Value(key3),Value(key4)),Value(key5),Value(key6)))
    {}

    template<class T_ARRAY>
    HASH(const T_ARRAY& array)
    {
        typedef typename T_ARRAY::ELEMENT T_ELEMENT;
        int i=1;int hash=array.Size()%2==0?missing_element_hash:Value(array(i++));
        for(;i<=array.Size()-1;i+=2) hash=triple_int_hash(hash,Value(array(i)),Value(array(i+1)));
        value=array.Size()==1?int_hash(hash):hash;
    }

private:
    static int Value(const int key)
    {return key;}

    static int Value(const HASH key)
    {return key.value;}

    template<class T>
    static int Value(const T& key)
    {return Value(Hash_Reduce(key));}

    unsigned int int_hash(unsigned int key) 
    {STATIC_ASSERT(sizeof(int)==4);
    key += ~(key << 15);
    key ^=  (key >> 10);
    key +=  (key << 3);
    key ^=  (key >> 6);
    key += ~(key << 11);
    key ^=  (key >> 16);
    return key;}

    unsigned int triple_int_hash(unsigned int a,unsigned int b,unsigned int c)
    {STATIC_ASSERT(sizeof(int)==4);
    a-=b;a-=c;a^=(c>>13);
    b-=c;b-=a;b^=(a<<8);
    c-=a;c-=b;c^=(b>>13);
    a-=b;a-=c;a^=(c>>12);
    b-=c;b-=a;b^=(a<<16);
    c-=a;c-=b;c^=(b>>5);
    a-=b;a-=c;a^=(c>>3);
    b-=c;b-=a;b^=(a<<10);
    c-=a;c-=b;c^=(b>>15);
    return c;}
};
//#####################################################################
// Function Hash_Reduce
//#####################################################################
inline int Hash_Reduce(const bool key){return key;}
inline int Hash_Reduce(const char key){return key;}
inline int Hash_Reduce(const unsigned char key){return key;}
inline int Hash_Reduce(const short key){return key;}
inline int Hash_Reduce(const unsigned short key){return key;}
inline int Hash_Reduce(const int key){return key;}
inline int Hash_Reduce(const unsigned int key){return key;}

inline int Hash_Reduce(const float key)
{STATIC_ASSERT(sizeof(float)==sizeof(int));
union {float f;int i;} raw;raw.f=key;return raw.i;}

inline HASH Hash_Reduce(const double key)
{STATIC_ASSERT(sizeof(double)==2*sizeof(int));
union {double d;int i[2];} raw;raw.d=key;
return HASH(1278312,raw.i[0],raw.i[1]);}

template<class T> inline HASH Hash_Reduce(const VECTOR<T,0>& key){return HASH();}
template<class T> inline HASH Hash_Reduce(const VECTOR<T,1>& key){return HASH(Hash_Reduce(key.x));}
template<class T> inline HASH Hash_Reduce(const VECTOR<T,2>& key){return HASH(key.x,key.y);}
template<class T> inline HASH Hash_Reduce(const VECTOR<T,3>& key){return HASH(key.x,key.y,key.z);}
template<class T> inline HASH Hash_Reduce(const VECTOR<T,4>& key){return HASH(key[1],key[2],key[3],key[4]);}
template<class T1,class T2> inline HASH Hash_Reduce(const PAIR<T1,T2>& key){return HASH(key.x,key.y);}
template<class T1,class T2,class T3> inline HASH Hash_Reduce(const TRIPLE<T1,T2,T3>& key){return HASH(key.x,key.y,key.z);}

template<class T> inline HASH Hash_Reduce(const COMPLEX<T>& key){return HASH(key.re,key.im);}
template<class T> inline HASH Hash_Reduce(const QUATERNION<T>& key){return HASH(key.s,key.v);}
template<class T> inline HASH Hash_Reduce(const ROTATION<VECTOR<T,1> >& key){return HASH();}
template<class T> inline HASH Hash_Reduce(const ROTATION<VECTOR<T,2> >& key){return Hash_Reduce(key.Complex());}
template<class T> inline HASH Hash_Reduce(const ROTATION<VECTOR<T,3> >& key){return Hash_Reduce(key.Quaternion());}

template<class TV> inline HASH Hash_Reduce(const RANGE<TV>& key){return HASH(key.min_corner,key.max_corner);}
template<class TV> inline HASH Hash_Reduce(const FRAME<TV>& key){return HASH(key.t,key.r);}
template<class TV> inline HASH Hash_Reduce(const TWIST<TV>& key){return HASH(key.linear,key.angular);}
template<int d> inline HASH Hash_Reduce(const FACE_INDEX<d>& key){return HASH(key.axis,key.index);}
template<int d> inline HASH Hash_Reduce(const SIDED_FACE_INDEX<d>& key){return HASH(key.side,key.axis,key.index);}

inline HASH Hash_Reduce(const char* key)
{return HASH(ARRAY_VIEW<const char>((int)strlen(key),key));}

inline HASH Hash_Reduce(const std::string& key)
{return HASH(ARRAY_VIEW<const char>((int)key.length(),key.c_str()));}

template<int s> inline HASH Hash_Reduce_Helper(const void* key);
template<> inline HASH Hash_Reduce_Helper<1>(const void* key){union {const void* p;int i;} raw;raw.p=key;return HASH(raw.i);}
template<> inline HASH Hash_Reduce_Helper<2>(const void* key){union {const void* p;int i[2];} raw;raw.p=key;return HASH(raw.i[0],raw.i[1]);}
inline HASH Hash_Reduce(const void* key){return Hash_Reduce_Helper<sizeof(void*)/sizeof(int)>(key);}

template<class T,class T_ARRAY> inline HASH Hash_Reduce(const ARRAY_BASE<T,T_ARRAY>& key)
{return HASH(key.Derived());}

template<class T,class T_ARRAY,int d> inline HASH Hash_Reduce(const ARRAY_BASE<T,T_ARRAY,VECTOR<int,d> >& key)
{return HASH(key.Derived().Domain_Indices(),key.array);}

template<class T,int length> inline HASH Hash_Reduce(const ARRAY<VECTOR<T,length>,FACE_INDEX<1> >& key){return Hash_Reduce(key.Component(1));}
template<class T,int length> inline HASH Hash_Reduce(const ARRAY<VECTOR<T,length>,FACE_INDEX<2> >& key){return HASH(key.Component(1),key.Component(2));}
template<class T,int length> inline HASH Hash_Reduce(const ARRAY<VECTOR<T,length>,FACE_INDEX<3> >& key){return HASH(key.Component(1),key.Component(2),key.Component(3));}

template<class ID,class T,int flags> inline HASH Hash_Reduce(const ELEMENT_ID<ID,T,flags>& id){return HASH(id.Value());}

//#####################################################################
// Function Hash
//#####################################################################
template<class T>
inline int Hash(const T& key)
{
    return HASH(Hash_Reduce(key)).value;
}
//#####################################################################
}
#endif
