//#########################################################################################################################################
// Copyright 2004-2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#########################################################################################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_KEY.h>
using namespace PhysBAM;
int OPENGL_KEY::MAX_KEY_INDEX=300;
int OPENGL_KEY::MAX_MODIFIER_INDEX=3;
//#####################################################################
// Function Name
//#####################################################################
std::string OPENGL_KEY::Name() const
{
    static const char* special_names[]={"F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12","Left","Right","Down","Up","PgDn","PgUp","Home","End","Ins","?"};
    std::string name;
    if(modifiers&ALT)  name+="Alt-";
    if(modifiers&CTRL) name+="^";
    if(key>=256) name+=std::string(special_names[(int)key-256]);
    else if (isprint(key)) name.append(1,(char)key);
    else if (key=='\t') name+="Tab";
    else name+="?";
    return name;
}
//#####################################################################
// Function Index
//#####################################################################
VECTOR<int,2> OPENGL_KEY::Index() const
{
    return VECTOR<int,2>(key,(int)modifiers);
}
//#####################################################################
// Function From_Glut_Key
//#####################################################################
OPENGL_KEY OPENGL_KEY::From_Glut_Key(unsigned char key,bool ctrl_pressed,bool alt_pressed)
{
    // need to correct for the fact that pressing ctrl alters the keycode for regular letters
    if('a'<=(key+'a'-1) && (key+'a'-1)<='z'){key+='a'-1; ctrl_pressed=true;}
    char modifiers=0x00;if(ctrl_pressed) modifiers|=CTRL;if(alt_pressed) modifiers|=ALT;
    return OPENGL_KEY(key,modifiers);
}
//#####################################################################
// Function From_Glut_Special_Key
//#####################################################################
OPENGL_KEY OPENGL_KEY::From_Glut_Special_Key(int key,bool ctrl_pressed,bool alt_pressed)
{
    char modifiers=0x00;if(ctrl_pressed) modifiers|=CTRL;if(alt_pressed) modifiers|=ALT;
    switch(key){
        case GLUT_KEY_F1: return OPENGL_KEY(F1,modifiers);
        case GLUT_KEY_F2: return OPENGL_KEY(F2,modifiers);
        case GLUT_KEY_F3: return OPENGL_KEY(F3,modifiers);
        case GLUT_KEY_F4: return OPENGL_KEY(F4,modifiers);
        case GLUT_KEY_F5: return OPENGL_KEY(F5,modifiers);
        case GLUT_KEY_F6: return OPENGL_KEY(F6,modifiers);
        case GLUT_KEY_F7: return OPENGL_KEY(F7,modifiers);
        case GLUT_KEY_F8: return OPENGL_KEY(F8,modifiers);
        case GLUT_KEY_F9: return OPENGL_KEY(F9,modifiers);
        case GLUT_KEY_F10: return OPENGL_KEY(F10,modifiers);
        case GLUT_KEY_F11: return OPENGL_KEY(F11,modifiers);
        case GLUT_KEY_F12: return OPENGL_KEY(F12,modifiers);
        case GLUT_KEY_LEFT: return OPENGL_KEY(LEFT,modifiers);
        case GLUT_KEY_UP: return OPENGL_KEY(UP,modifiers);
        case GLUT_KEY_RIGHT: return OPENGL_KEY(RIGHT,modifiers);
        case GLUT_KEY_DOWN: return OPENGL_KEY(DOWN,modifiers);
        case GLUT_KEY_PAGE_UP: return OPENGL_KEY(PAGE_UP,modifiers);
        case GLUT_KEY_PAGE_DOWN: return OPENGL_KEY(PAGE_DOWN,modifiers);
        case GLUT_KEY_HOME: return OPENGL_KEY(HOME,modifiers);
        case GLUT_KEY_END: return OPENGL_KEY(END,modifiers);
        case GLUT_KEY_INSERT: return OPENGL_KEY(INSERT,modifiers);
        default: return OPENGL_KEY();}
}
//#####################################################################
// Function From_String
//#####################################################################
OPENGL_KEY OPENGL_KEY::From_String(const std::string& key_string)
{
    unsigned int dummy_index=0;return From_String(key_string,dummy_index);
}
//#####################################################################
// Function From_String
//#####################################################################
// i is the index in key_string at which to start parsing, and is updated to be the index pointing to the start of the next token to parse
// returned key OPENGL_KEY(UNKNOWN) indicates an error
OPENGL_KEY OPENGL_KEY::From_String(const std::string& key_string,unsigned int& i)
{
    static const char* key_names[]={"F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12","LEFT","RIGHT","DOWN","UP","PAGE_DOWN","PAGE_UP","HOME","END","INSERT",0};

    if(i>=key_string.length()) return OPENGL_KEY(UNKNOWN);
    else if(key_string[i]=='\\'){ // suppress special meaning of what follows
        if(i+1<key_string.length()){i++;return OPENGL_KEY(key_string[i++]);} else return OPENGL_KEY(UNKNOWN);}
    else if(key_string[i]=='^'){ // ctrl-key sequence
        if(i+1<key_string.length()){i++;return OPENGL_KEY(key_string[i++],OPENGL_KEY::CTRL);} else return OPENGL_KEY(UNKNOWN);}
    else if(key_string[i]=='<'){ // special named key sequence
        std::string::size_type end_index=key_string.find('>',i+1);
        if(end_index==std::string::npos){i=(unsigned int)key_string.size();return OPENGL_KEY(UNKNOWN);}
        else{std::string key_name=STRING_UTILITIES::toupper(key_string.substr(i+1,end_index-i-1));i=(unsigned int)end_index+1;
            for(int j=0;key_names[j];j++) if(key_name==std::string(key_names[j])) return OPENGL_KEY((OPENGL_KEY::SPECIAL_KEY)((int)OPENGL_KEY::F1+j));
            return OPENGL_KEY(UNKNOWN);}}
    else{return key_string[i++];}
}
//#####################################################################
// Function Parse_Key_Sequence
//#####################################################################
void OPENGL_KEY::Parse_Key_Sequence(const std::string& key_string,ARRAY<OPENGL_KEY>& key_list)
{
    key_list.Remove_All();
    unsigned int i=0;
    while(i<key_string.length()){
        key_list.Append(From_String(key_string,i));
        if(key_list(key_list.m).key==OPENGL_KEY::UNKNOWN){LOG::cerr << "Unknown key sequence reached at index " << i << " of '" << key_string << "'" << std::endl;break;}}
}
//#####################################################################
