//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEBUG_UTILITIES
//#####################################################################
#ifndef __DEBUG_UTILITIES__
#define __DEBUG_UTILITIES__

#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <string>
#include <typeinfo>

#define PHYSBAM_DEBUG_FUNCTION_NAME ((const char*)__FUNCTION__) // cast to const char* to work around error in noreturn

#define PHYSBAM_WARN_IF_NOT_OVERRIDDEN() \
    do{static bool __first_time__=true;if(__first_time__){PhysBAM::DEBUG_UTILITIES::Warn_If_Not_Overridden(__FUNCTION__,__FILE__,__LINE__,typeid(*this));__first_time__=false;}}while(0)

#define PHYSBAM_WARNING(message) \
    do{static bool __first_time__=true;if(__first_time__){PhysBAM::DEBUG_UTILITIES::Warning((message),__FUNCTION__,__FILE__,__LINE__);__first_time__=false;}}while(0)

#define PHYSBAM_FUNCTION_IS_NOT_DEFINED() \
    PhysBAM::DEBUG_UTILITIES::Function_Is_Not_Defined(PHYSBAM_DEBUG_FUNCTION_NAME,__FILE__,__LINE__,typeid(*this))

#define PHYSBAM_NOT_IMPLEMENTED(...) \
    PhysBAM::DEBUG_UTILITIES::Not_Implemented(PHYSBAM_DEBUG_FUNCTION_NAME,__FILE__,__LINE__,##__VA_ARGS__)

#define PHYSBAM_FATAL_ERROR(...) \
    PhysBAM::DEBUG_UTILITIES::Fatal_Error(PHYSBAM_DEBUG_FUNCTION_NAME,__FILE__,__LINE__,##__VA_ARGS__)

#define PHYSBAM_ASSERT(condition,...) \
    if(condition){}else PhysBAM::DEBUG_UTILITIES::Assertion_Failed(PHYSBAM_DEBUG_FUNCTION_NAME,__FILE__,__LINE__,#condition,##__VA_ARGS__)

#ifdef NDEBUG
#   define PHYSBAM_DEBUG_ONLY(...)
#else
#   define PHYSBAM_DEBUG_ONLY(...) __VA_ARGS__
#endif

namespace PhysBAM{
namespace DEBUG_UTILITIES{
class DEBUG_GLOBAL_VARIABLES{
    public:
        static int int_variable;
};
void Debug_Breakpoint();
void Warn_If_Not_Overridden(const char* function,const char* file,unsigned int line,const std::type_info& type);
void Warning(const std::string& message,const char* function,const char* file,unsigned int line);
void PHYSBAM_NORETURN(Function_Is_Not_Defined(const char* function,const char* file,unsigned int line,const std::type_info& type));
void PHYSBAM_NORETURN(Not_Implemented(const char* function,const char* file,unsigned int line)); // three different versions to minimize caller code
void PHYSBAM_NORETURN(Not_Implemented(const char* function,const char* file,unsigned int line,const char* message));
void PHYSBAM_NORETURN(Not_Implemented(const char* function,const char* file,unsigned int line,const std::string& message));
void PHYSBAM_NORETURN(Fatal_Error(const char* function,const char* file,unsigned int line)); // three different versions to minimize caller code
void PHYSBAM_NORETURN(Fatal_Error(const char* function,const char* file,unsigned int line,const char* message));
void PHYSBAM_NORETURN(Fatal_Error(const char* function,const char* file,unsigned int line,const std::string& message));
void PHYSBAM_NORETURN(Assertion_Failed(const char* function,const char* file,unsigned int line,const char* condition)); // three different versions to minimize caller code
void PHYSBAM_NORETURN(Assertion_Failed(const char* function,const char* file,unsigned int line,const char* condition,const char* message));
void PHYSBAM_NORETURN(Assertion_Failed(const char* function,const char* file,unsigned int line,const char* condition,const std::string& message));
}
}
#endif
