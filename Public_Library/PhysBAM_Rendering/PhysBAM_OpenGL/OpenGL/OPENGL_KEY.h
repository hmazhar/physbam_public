//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_KEY
//#####################################################################
#ifndef __OPENGL_KEY__
#define __OPENGL_KEY__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h> // just so we get gl in the right order
#include <string>
namespace PhysBAM{

class OPENGL_KEY
{
public:
    static int MAX_KEY_INDEX,MAX_MODIFIER_INDEX;
    enum MODIFIER{CTRL=0x01,ALT=0x02}; // update MAX_MODIFIER_INDEX if changed!
    enum SPECIAL_KEY{F1=256,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,LEFT,RIGHT,DOWN,UP,PAGE_DOWN,PAGE_UP,HOME,END,INSERT,UNKNOWN}; // update MAX_KEY_INDEX if changed!
    int key; // 0-255 for normal keys, 256- for special keys (using enums)
    char modifiers; // a bit-OR of the MODIFIER enums
   
    OPENGL_KEY()
        :key((int)UNKNOWN),modifiers(0x00)
    {}

    OPENGL_KEY(const unsigned char key_input,const char modifiers_input=0x00)
        :key((int)key_input),modifiers(modifiers_input)
    {}

    OPENGL_KEY(const SPECIAL_KEY special_key,const char modifiers_input=0x00)
        :key((int)special_key),modifiers(modifiers_input)
    {}
    
    bool operator==(const OPENGL_KEY& other_key)
    {return key==other_key.key && modifiers==other_key.modifiers;}

//#########################################################################################################################################
    std::string Name() const;
    VECTOR<int,2> Index() const; // Each key gets a unique 2d index for fast lookups
    static OPENGL_KEY From_Glut_Key(unsigned char key,bool ctrl_pressed=false,bool alt_pressed=false);
    static OPENGL_KEY From_Glut_Special_Key(int key,bool ctrl_pressed=false,bool alt_pressed=false);
    static OPENGL_KEY From_String(const std::string& key_string);
    static OPENGL_KEY From_String(const std::string& key_string,unsigned int& i);
    static void Parse_Key_Sequence(const std::string& key_string,ARRAY<OPENGL_KEY>& key_list);
//#########################################################################################################################################
};
}
#endif
