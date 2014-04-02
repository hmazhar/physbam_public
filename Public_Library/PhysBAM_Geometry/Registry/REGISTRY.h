//#####################################################################
// Copyright 2006, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACTORY
//#####################################################################
#ifndef __REGISTRY__
#define __REGISTRY__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#endif
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <string>
namespace PhysBAM{

template<class T_BASE_OBJECT>
struct FACTORY_BASE
{
    typedef typename T_BASE_OBJECT::VECTOR_T TV;typedef typename TV::SCALAR T;

    virtual ~FACTORY_BASE(){}
    virtual T_BASE_OBJECT* Create() const=0;
    virtual T_BASE_OBJECT* Create(GEOMETRY_PARTICLES<TV>& particles) const=0;
};

template<class T_OBJECT,class T_BASE_OBJECT>
struct FACTORY:public FACTORY_BASE<T_BASE_OBJECT>
{
    typedef typename T_BASE_OBJECT::VECTOR_T TV;typedef typename TV::SCALAR T;

    typedef T_OBJECT* CREATE_RAW();
    typedef T_OBJECT* CREATE_WITH_PARTICLES(GEOMETRY_PARTICLES<TV>&);

    template<class T_CREATE> T_CREATE*
    Get_Create_Helper(...) const
    {return 0;}

    template<class T_CREATE> T_CREATE*
    Get_Create_Helper(T_CREATE* create) const
    {return create;}

    template<class T_CREATE> T_CREATE*
    Get_Create() const
    {return Get_Create_Helper<T_CREATE>(&T_OBJECT::Create);}

    T_OBJECT* Create() const PHYSBAM_OVERRIDE
    {CREATE_RAW* create=Get_Create<CREATE_RAW>();
    return create?create():0;}

    T_OBJECT* Create(GEOMETRY_PARTICLES<TV>& particles) const PHYSBAM_OVERRIDE
    {CREATE_WITH_PARTICLES* create=Get_Create<CREATE_WITH_PARTICLES>();
    return create?create(particles):0;}
};

//#####################################################################
// Class REGISTRY
//#####################################################################
template<class T_BASE_OBJECT,class T_NAME,class T_DERIVED_REGISTRY>
class REGISTRY
{
private:
    HASHTABLE<T_NAME,FACTORY_BASE<T_BASE_OBJECT>*> name_registry;
    HASHTABLE<std::string,FACTORY_BASE<T_BASE_OBJECT>*> extension_registry;

    static T_DERIVED_REGISTRY& Singleton()
    {static T_DERIVED_REGISTRY registry;return registry;}

protected:
    REGISTRY()
    {}

    ~REGISTRY()
    {name_registry.Delete_Pointers_Stored_In_Table();}

    template<class T_OBJECT> static T_OBJECT* Create_Representative()
    {return T_OBJECT::Create();}

public:
    // Register cannot be called from the constructor of a derived class, since this would cause a recursive call to Singleton
    template<class T_OBJECT> static void Register()
    {REGISTRY& registry=Singleton();
    T_NAME name=T_OBJECT::Static_Name();std::string extension=T_OBJECT::Static_Extension();
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    //LOG::cout<<"Registering "<<name<<" to extension "<<extension<<std::endl;
#endif

    FACTORY_BASE<T_BASE_OBJECT>* factory=new FACTORY<T_OBJECT,T_BASE_OBJECT>();
    if(registry.name_registry.Contains(name)) PHYSBAM_FATAL_ERROR();
    registry.name_registry.Insert(name,factory);
    if(!extension.empty()){    
        if(registry.extension_registry.Contains(extension)) PHYSBAM_FATAL_ERROR();
        registry.extension_registry.Insert(extension,factory);}}

    static FACTORY_BASE<T_BASE_OBJECT>* Name_To_Factory(const T_NAME& name)
    {REGISTRY& registry=Singleton();
    FACTORY_BASE<T_BASE_OBJECT>* factory;
    if(!registry.name_registry.Get(name,factory)){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cerr<<name<<" is not a registered object name.\nRegisted names are "<<registry.name_registry<<std::endl;
#endif
        PHYSBAM_FATAL_ERROR();}
    return factory;}

    static FACTORY_BASE<T_BASE_OBJECT>* Extension_To_Factory(const std::string& extension)
    {REGISTRY& registry=Singleton();
    FACTORY_BASE<T_BASE_OBJECT>* factory;
    if(!registry.extension_registry.Get(extension,factory)){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cerr<<extension<<" is not a registered object extension.\nRegisted extensions are "<<registry.extension_registry<<std::endl;
#endif
        PHYSBAM_FATAL_ERROR();}
    return factory;}

//#####################################################################
};
}
#endif
