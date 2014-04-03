//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEBUG_UTILITIES
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Utilities/EXCEPTIONS.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <cassert>
#include <stdexcept>
#include <typeinfo>
#if defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__)
#include <csignal>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace PhysBAM{
namespace DEBUG_UTILITIES{
int DEBUG_GLOBAL_VARIABLES::int_variable=0;
#if defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__)
//#####################################################################
// Function Debug_Breakpoint Linux
//#####################################################################
void Debug_Breakpoint()
{
    kill(getpid(),SIGINT); // if you use this you need to step out of the signal handler to get a non-corrupt stack
}
#else
//#####################################################################
// Function Debug_Breakpoint Default
//#####################################################################
void Debug_Breakpoint()
{
    assert(false);
}
#endif
//#####################################################################
// Function Warn_If_Not_Overridden
//#####################################################################
void Warn_If_Not_Overridden(const char* function,const char* file,unsigned int line,const std::type_info& type)
{
    LOG::cerr<<STRING_UTILITIES::string_sprintf("*** PHYSBAM_WARNING: %s:%s:%d: Function not overridden by %s",file,function,line,type.name())<<std::endl;
}
//#####################################################################
// Function Warning
//#####################################################################
void Warning(const std::string& message,const char* function,const char* file,unsigned int line)
{
    LOG::cerr<<STRING_UTILITIES::string_sprintf("*** PHYSBAM_WARNING: %s:%s:%d: %s",file,function,line,message.c_str())<<std::endl;
}
//#####################################################################
// Function Function_Is_Not_Defined
//#####################################################################
void Function_Is_Not_Defined(const char* function,const char* file,unsigned int line,const std::type_info& type)
{
    std::string error=STRING_UTILITIES::string_sprintf("%s:%s:%d: Function not defined by %s",file,function,line,type.name());
    LOG::cout<<std::flush;LOG::cerr<<"\n";
    PROCESS_UTILITIES::Backtrace();
    LOG::cerr<<"\n*** ERROR: "<<error<<'\n'<<std::endl;
    throw std::runtime_error(error);
}
//#####################################################################
// Function Not_Implemented
//#####################################################################
void Not_Implemented(const char* function,const char* file,unsigned int line)
{
    Not_Implemented(function,file,line,"Something");
}
void Not_Implemented(const char* function,const char* file,unsigned int line,const char* message)
{
    Not_Implemented(function,file,line,std::string(message));
}
void Not_Implemented(const char* function,const char* file,unsigned int line,const std::string& message)
{
    std::string error=STRING_UTILITIES::string_sprintf("%s:%s:%d: Not implemented: %s",file,function,line,message.c_str());
    LOG::cout<<std::flush;LOG::cerr<<"\n";
    PROCESS_UTILITIES::Backtrace();
    LOG::cerr<<"\n*** ERROR: "<<error<<'\n'<<std::endl;
    PROCESS_UTILITIES::Set_Backtrace(false);
    throw NOT_IMPLEMENTED_ERROR(error);
}
//#####################################################################
// Function Fatal_Error
//#####################################################################
void Fatal_Error(const char* function,const char* file,unsigned int line)
{
    Fatal_Error(function,file,line,"Fatal error");
}
void Fatal_Error(const char* function,const char* file,unsigned int line,const char* message)
{
    Fatal_Error(function,file,line,std::string(message));
}
void Fatal_Error(const char* function,const char* file,unsigned int line,const std::string& message)
{
    std::string error=STRING_UTILITIES::string_sprintf("%s:%s:%d: %s",file,function,line,message.c_str());
    LOG::cout<<std::flush;LOG::cerr<<"\n";
    PROCESS_UTILITIES::Backtrace();
    LOG::cerr<<"\n*** ERROR: "<<error<<'\n'<<std::endl;
    PROCESS_UTILITIES::Set_Backtrace(false);
    throw std::runtime_error(error);
}
//#####################################################################
// Function Assertion_Failed
//#####################################################################
void Assertion_Failed(const char* function,const char* file,unsigned int line,const char* condition)
{
    Assertion_Failed(function,file,line,condition,"Assertion failed");
}
void Assertion_Failed(const char* function,const char* file,unsigned int line,const char* condition,const char* message)
{
    Assertion_Failed(function,file,line,condition,std::string(message));
}
void Assertion_Failed(const char* function,const char* file,unsigned int line,const char* condition,const std::string& message)
{
    std::string error=STRING_UTILITIES::string_sprintf("%s:%s:%d: %s, condition = %s",file,function,line,message.c_str(),condition);
    LOG::cout<<std::flush;LOG::cerr<<"\n";
    PROCESS_UTILITIES::Backtrace();
    LOG::cerr<<"\n*** ERROR: "<<error<<'\n'<<std::endl;
    PROCESS_UTILITIES::Set_Backtrace(false);
    throw ASSERTION_ERROR(error);
}
//#####################################################################
}
}
