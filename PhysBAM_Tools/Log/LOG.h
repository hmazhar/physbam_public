//#####################################################################
// Copyright 2004-2008, Geoffrey Irving, Frank Losasso, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOG
//##################################################################### 
#ifndef __LOG__
#define __LOG__

#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <cassert>
#include <ostream>
#include <sstream>
namespace PhysBAM{

#define LOG_REAL LOG

namespace LOG_NULL{
struct log_null_class {
    template<class T> log_null_class& operator<<(const T&){return *this;}
    log_null_class& operator<<(std::ostream& (*)(std::ostream&)){return *this;}
    template<class T> void flags(const T&);
    template<class T> void width(const T&);
    template<class A,class B> void Copy_Log_To_File(const A&,const B&){}
};
extern log_null_class cout;
extern log_null_class cerr;
template<class T> inline void Time(const T&){}
template<class A,class B> inline void Stat(const A&,const B&){}
inline void Stop_Time(){}
inline void Finish_Logging(){}
template<class A,class B,class C,class D> inline void Initialize_Logging(const A&,const B&,const C&,const D&){}
inline log_null_class* Instance(){return &cout;}

struct SCOPE
{
    template<class A> SCOPE(const A&){}
    template<class A,class B> SCOPE(const A&,const B&){}
    template<class A,class B,class C> SCOPE(const A&,const B&,const C&){}
    template<class A,class B,class C,class D> SCOPE(const A&,const B&,const C&,const D&){}
    template<class A,class B,class C,class D,class E> SCOPE(const A&,const B&,const C&,const D&,const E&){}
    void Pop(){}
};
}

namespace LOG_REAL{

class LOG_ENTRY;
class LOG_SCOPE;

class LOG_CLASS
{
    friend class LOG_ENTRY;
    friend class LOG_SCOPE;
    friend class LOG_COUT_BUFFER;
    friend class LOG_CERR_BUFFER;
    friend void Reset();
    friend void Dump_Log();

    TIMER* timer_singleton;
    int timer_id;
    bool suppress_cout;
    bool suppress_cerr;
public:
    bool suppress_timing;
    FILE* log_file;
    int verbosity_level;
    bool log_file_temporary;
    bool xml;

    LOG_ENTRY* root;
    LOG_ENTRY* current_entry;

    LOG_CLASS(const bool suppress_cout,const bool suppress_cerr,const bool suppress_timing,const int verbosity_level,const bool cache_initial_output);
    ~LOG_CLASS();

    static void Push_Scope(const std::string& scope_identifier,const std::string& scope_name);
    static void Pop_Scope();
public:
    static void Time_Helper(const std::string& label);
    void Copy_Log_To_File(const std::string& filename,const bool append);

//##################################################################### 
};

class SCOPE:private NONCOPYABLE
{
    bool active;
public:
    SCOPE()
        :active(false)
    {}

    SCOPE(const std::string& scope_identifier)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,scope_identifier);
    }

    SCOPE(const std::string& scope_identifier,const std::string& scope_name)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,scope_name);
    }

    template<class T1>
    SCOPE(const std::string& scope_identifier,const std::string& format,const T1& d1)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,STRING_UTILITIES::string_sprintf(format.c_str(),d1));
    }

    template<class T1,class T2>
    SCOPE(const std::string& scope_identifier,const std::string& format,const T1& d1,const T2& d2)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,STRING_UTILITIES::string_sprintf(format.c_str(),d1,d2));
    }

    template<class T1,class T2,class T3>
    SCOPE(const std::string& scope_identifier,const std::string& format,const T1& d1,const T2& d2,const T3& d3)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,STRING_UTILITIES::string_sprintf(format.c_str(),d1,d2,d3));
    }

    template<class T1,class T2,class T3,class T4>
    SCOPE(const std::string& scope_identifier,const std::string& format,const T1& d1,const T2& d2,const T3& d3,const T4& d4)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,STRING_UTILITIES::string_sprintf(format.c_str(),d1,d2,d3,d4));
    }

    template<class T1,class T2,class T3,class T4,class T5>
    SCOPE(const std::string& scope_identifier,const std::string& format,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5)
        :active(true)
    {
        LOG_CLASS::Push_Scope(scope_identifier,STRING_UTILITIES::string_sprintf(format.c_str(),d1,d2,d3,d4,d5));
    }

    ~SCOPE();
    
    void Push(const std::string& scope_identifier,const std::string& scope_name)
    {assert(!active);active=true;LOG_CLASS::Push_Scope(scope_identifier,scope_name);}

    template<class T1>
    void Push(const std::string& scope_identifier,const std::string& format,const T1& d1)
    {Push(scope_identifier,STRING_UTILITIES::string_sprintf(format,d1));}

    void Pop()
    {assert(active);active=false;LOG_CLASS::Pop_Scope();}
//##################################################################### 
};

// These next few lines are important to ensure no static data from LOG.cpp is accessed for DLLs
LOG_CLASS* Instance();
std::ostream& cout_Helper();
std::ostream& cerr_Helper();
namespace{
    static std::ostream& cout PHYSBAM_UNUSED =::PhysBAM::LOG_REAL::cout_Helper();
    static std::ostream& cerr PHYSBAM_UNUSED =::PhysBAM::LOG_REAL::cerr_Helper();
}

void Initialize_Logging(const bool suppress_cout_input=false,const bool suppress_timing_input=false,const int verbosity_level_input=1<<30,const bool cache_initial_output=false);
void Finish_Logging();
void Stop_Time();

//#####################################################################
// Function Stat
//#####################################################################
void Stat_Helper(const std::string& label,const std::stringstream& s);
template<class T_VALUE> void
Stat(const std::string& label,const T_VALUE& value)
{std::stringstream s;s<<value;Stat_Helper(label,s);}
void Reset(); 
void Dump_Log();

inline void Time(const std::string& format)
{if(Instance()->suppress_timing) return;
LOG_CLASS::Time_Helper(format);}

//##################################################################### 
}
}
#endif
