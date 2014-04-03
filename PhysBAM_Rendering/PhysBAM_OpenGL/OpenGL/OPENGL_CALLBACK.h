//#####################################################################
// Copyright 2002-2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_CALLBACK
//##################################################################### 
#ifndef __OPENGL_CALLBACK__
#define __OPENGL_CALLBACK__
#include <iostream>
#include <string>
// A convenience macro
//
// In a class you put e.g.:
//
// class FOO
// {
//      void Toggle_Draw();
//      DEFINE_CALLBACK_CREATOR(FOO, Toggle_Draw);
// };
//
// And you define the Toggle_Draw function.
// Then to bind this to a key (assuming you have FOO *foo) you use
//
// opengl_world.Bind('d', foo->Toggle_Draw_CB());
//
#define DEFINE_CALLBACK_CREATOR(classname, func) \
    class classname##func : public PhysBAM::OPENGL_CALLBACK { \
        public: \
            classname##func(classname *obj, const std::string &help_string) : obj(obj), help_string(help_string) {} \
            void operator() () { obj->func(); } \
            void Print(std::ostream &out) { out << help_string; } \
        private: \
             classname *obj; \
             std::string help_string; }; \
    PhysBAM::OPENGL_CALLBACK *func##_CB(const std::string &help_string = "Unknown") \
    { return new classname##func(this, help_string); }

namespace PhysBAM{

class OPENGL_CALLBACK
{
public:
    virtual ~OPENGL_CALLBACK();
    virtual void operator() () = 0;
    virtual void Print(std::ostream& out);
//#####################################################################
};
}
#endif

