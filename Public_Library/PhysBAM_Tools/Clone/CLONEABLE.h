//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CLONEABLE 
//#####################################################################
#ifndef __CLONEABLE__
#define __CLONEABLE__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <typeinfo>
namespace PhysBAM{

class CLONEABLE_BASE; // cloneable abstract base classes inherit from this directly
template<class T> class CLONE_ARRAY;

template<class T> struct IS_CLONEABLE:public AND<
    IS_BASE_OF<CLONEABLE_BASE,T>::value, // must eventually derive from CLONEABLE_BASE
    IS_CONVERTIBLE<T*,CLONEABLE_BASE*>::value>{}; // ensure derivation isn't ambiguous

//#####################################################################
// Class CLONEABLE_BASE
//#####################################################################
class CLONEABLE_BASE
{
public:
    virtual ~CLONEABLE_BASE()
    {}

    CLONEABLE_BASE* Clone() const
    {return Clone_Implementation();}

    CLONEABLE_BASE* Clone_Default() const // creates a default constructed copy of the same type
    {return Clone_Default_Implementation();}

protected:
    template<class T> friend class CLONE_ARRAY;

    virtual CLONEABLE_BASE* Clone_Implementation() const=0;
    virtual CLONEABLE_BASE* Clone_Default_Implementation() const=0;
    virtual size_t Sizeof_Clone() const=0;
    virtual CLONEABLE_BASE* Placement_Clone(void* memory) const=0;
};
//#####################################################################
// Class CLONEABLE_ABSTRACT
//#####################################################################
template<class T_DERIVED,class T_BASE=CLONEABLE_BASE>
class CLONEABLE_ABSTRACT:public T_BASE
{
    STATIC_ASSERT(IS_CLONEABLE<T_BASE>::value);
    using T_BASE::Clone_Implementation;using T_BASE::Clone_Default_Implementation;
    template<class T> friend class CLONE_ARRAY;
public:

    T_DERIVED* Clone() const
    {return dynamic_cast<T_DERIVED*>(Clone_Implementation());}

    T_DERIVED* Clone_Default() const // creates a default constructed copy of the same type
    {return dynamic_cast<T_DERIVED*>(Clone_Default_Implementation());}
};
//#####################################################################
// Class CLONEABLE
//#####################################################################
template<class T_DERIVED,class T_BASE=CLONEABLE_BASE>
class CLONEABLE:public T_BASE
{
    STATIC_ASSERT(IS_CLONEABLE<T_BASE>::value);
public:

    T_DERIVED* Clone() const
    {return dynamic_cast<T_DERIVED*>(Clone_Implementation());}

    T_DERIVED* Clone_Default() const // creates a default constructed copy of the same type
    {return dynamic_cast<T_DERIVED*>(Clone_Default_Implementation());}

private:
    template<class T> friend class CLONE_ARRAY;

    virtual CLONEABLE_BASE* Clone_Implementation() const
    {PHYSBAM_ASSERT(typeid(*this)==typeid(T_DERIVED)); // avoid slicing errors
    T_DERIVED* clone=new T_DERIVED();
    clone->Clone_Helper(dynamic_cast<const T_DERIVED&>(*this));
    return clone;}

    virtual CLONEABLE_BASE* Clone_Default_Implementation() const
    {PHYSBAM_ASSERT(typeid(*this)==typeid(T_DERIVED)); // avoid slicing errors
    return new T_DERIVED();}

    virtual size_t Sizeof_Clone() const
    {PHYSBAM_ASSERT(typeid(*this)==typeid(T_DERIVED)); // avoid horrible memory corruption errors
    return sizeof(T_DERIVED);}

    virtual CLONEABLE_BASE* Placement_Clone(void* memory) const
    {PHYSBAM_ASSERT(typeid(*this)==typeid(T_DERIVED)); // avoid slicing errors
    T_DERIVED* clone=new(memory) T_DERIVED();
    clone->Clone_Helper(dynamic_cast<const T_DERIVED&>(*this));
    return clone;}
};
//#####################################################################
}
#endif
