//#####################################################################
// Copyright 2002-2010, Eran Guendelman, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_BASIC_CALLBACKS
//#####################################################################
#ifndef __OPENGL_BASIC_CALLBACKS__
#define __OPENGL_BASIC_CALLBACKS__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
namespace PhysBAM{

//#####################################################################
// Class OPENGL_CALLBACK_FULLSCREEN
//#####################################################################
// a helper callback class that is glut specific
class OPENGL_CALLBACK_FULLSCREEN:public OPENGL_CALLBACK
{
    bool fullscreen_mode;
    int *width,*height,saved_width,saved_height;
public:
    OPENGL_CALLBACK_FULLSCREEN(int *width_input,int *height_input)
        :fullscreen_mode(false),width(width_input),height(height_input),saved_width(*width_input),saved_height(*height_input)
    {}

    void operator()()
    {if(fullscreen_mode){fullscreen_mode=false;*width=saved_width;*height=saved_height;glutReshapeWindow(saved_width,saved_height);}
    else{fullscreen_mode=true;saved_width=*width;saved_height=*height;glutFullScreen();}}

    void Print(std::ostream& out)
    {out<<"Toggle fullscreen";}
};
//#####################################################################
// Class OPENGL_CALLBACK_TOGGLE
//#####################################################################
class OPENGL_CALLBACK_TOGGLE:public OPENGL_CALLBACK
{
    bool *state;
    std::string help_string;
public:
    OPENGL_CALLBACK_TOGGLE(bool *state_input, const std::string &help_string_input="Toggle state") 
        : state(state_input), help_string(help_string_input) {}
    void operator() () PHYSBAM_OVERRIDE {*state=!*state;}
    void Print(std::ostream& out) { out << help_string << " (" << (*state?"1":"0") << ")"; }
};
//#####################################################################
// Class OPENGL_CALLBACK_CYCLE
//#####################################################################
class OPENGL_CALLBACK_CYCLE:public OPENGL_CALLBACK
{
    int *state;
    int min_value, max_value;
    std::string help_string;
public:
    OPENGL_CALLBACK_CYCLE(int *state, int min_value, int max_value, const std::string &help_string_input="Cycle state") 
        : state(state), min_value(min_value), max_value(max_value), help_string(help_string_input) {}
    void operator() () PHYSBAM_OVERRIDE { if (++(*state) > max_value) *state = min_value; }
    void Print(std::ostream& out) { out << help_string << " (" << *state << ")"; }
};
//#####################################################################
// Class OPENGL_CALLBACK_SET_VALUE
//#####################################################################
template<class T>
class OPENGL_CALLBACK_SET_VALUE:public OPENGL_CALLBACK
{
    T *state, value;
    std::string help_string;
public:
    OPENGL_CALLBACK_SET_VALUE(T *state_input, T value_input, const std::string &help_string_input="Set value") 
        : state(state_input),value(value_input),help_string(help_string_input) {}
    void operator() () PHYSBAM_OVERRIDE {*state=value;}
    void Print(std::ostream& out) { out << help_string << " (" << *state << "=" << value << ")"; }
};
//#####################################################################
// Class OPENGL_CALLBACK_ZOOM
//#####################################################################
// actually a misnomer - this is panning in/out, not zooming in the true sense of the word
class OPENGL_CALLBACK_ZOOM:public OPENGL_CALLBACK
{
    const float factor;
    float *camera_distance, *nearclip, *farclip;
public:
    OPENGL_CALLBACK_ZOOM(const float factor_input,float *camera_distance_input,float *nearclip_input,float *farclip_input)
        : factor(factor_input),camera_distance(camera_distance_input),nearclip(nearclip_input),farclip(farclip_input)
    {}

    void operator() () PHYSBAM_OVERRIDE {*camera_distance*=factor;*nearclip*=factor;*farclip*=factor;}
    void Print(std::ostream& out) { out << "Zoom by " << factor; }
};

//=====================================================================
// OPENGL_WORLD specific callbacks
//=====================================================================

//#####################################################################
// Class OPENGL_CALLBACK_SAVE_SCREEN
//#####################################################################
class OPENGL_CALLBACK_SAVE_SCREEN:public OPENGL_CALLBACK
{
    OPENGL_WORLD &world;
    const std::string basename;
    int counter;
public:
    OPENGL_CALLBACK_SAVE_SCREEN(OPENGL_WORLD &world_input,const std::string& basename_input="opengl",int counter_start=0)
        : world(world_input),basename(basename_input),counter(counter_start)
    {}

    void operator() () PHYSBAM_OVERRIDE
    {
        std::string defaultname=STRING_UTILITIES::string_sprintf("%s%03d.ppm",basename.c_str(),counter);
        LOG::cout<<"Filename for saving image in .ppm or .jpg (if supported) formats (return for default, - for cancel)\n";
        LOG::cout<<"[%s]: "<<defaultname;
        std::string filename;getline(std::cin,filename);
        size_t len=filename.length();if(filename[len-1]=='\n'){--len;filename=filename.substr(0,len);}
        if(len==0)filename=defaultname;
        else if(filename=="-") return;
        else{} // should parse user input to make a better default basename next time
        LOG::cout<<"saving to %s..."<<filename;
        world.Save_Screen(filename,true);
        LOG::cout<<"\n";
        ++counter;
    }
    void Print(std::ostream& out) { out << "Save screen"; }
};
//#####################################################################
// Class OPENGL_CALLBACK_SAVE_TO_EPS
//#####################################################################
class OPENGL_CALLBACK_SAVE_TO_EPS:public OPENGL_CALLBACK
{
    OPENGL_WORLD &world;
    const std::string basename;
    int counter;
public:
    OPENGL_CALLBACK_SAVE_TO_EPS(OPENGL_WORLD &world_input,const std::string& basename_input="opengl",int counter_start=0)
        : world(world_input),basename(basename_input),counter(counter_start)
    {}

    void operator() () PHYSBAM_OVERRIDE
    {
        std::string filename=STRING_UTILITIES::string_sprintf("%s%03d.eps",basename.c_str(),counter);
        opengl_eps_output = new OPENGL_EPS_OUTPUT<float>(filename);

        LOG::cout<<"saving to %s..."<<filename<<std::endl;

        world.Render_World(false);
        glFinish();
        delete opengl_eps_output;opengl_eps_output=0;

        ++counter;
    }
    void Print(std::ostream& out) { out << "Save screen to eps"; }
};
//#####################################################################
// Class OPENGL_CALLBACK_MOVE_TARGET
//#####################################################################
class OPENGL_CALLBACK_MOVE_TARGET:public OPENGL_CALLBACK
{
    OPENGL_WORLD &world;
    const VECTOR<float,3> motion_step;
    const float *motion_factor;
    VECTOR<float,3> *target;
public:
    OPENGL_CALLBACK_MOVE_TARGET(OPENGL_WORLD& world_input,const VECTOR<float,3>& motion_step_input,const float* motion_factor_input,VECTOR<float,3>* target_input)
        :world(world_input),motion_step(motion_step_input),motion_factor(motion_factor_input),target(target_input)
    {}

    void operator() ()  PHYSBAM_OVERRIDE
    {
        VECTOR<float,3> view_forward,view_up,view_right;
        world.Get_View_Frame(view_forward,view_up,view_right);
        (*target)+=*motion_factor*(motion_step.x*view_right+motion_step.y*view_up+motion_step.z*view_forward);
        world.Set_View_Target_Timer(1);
    }
    void Print(std::ostream& out) PHYSBAM_OVERRIDE
    {if(motion_step!=VECTOR<float,3>()) out<<"Move target "<<motion_step;else out<<"View target";}
};
//#####################################################################
// Class OPENGL_CALLBACK_TOGGLE_SMOOTH_SHADING
//#####################################################################
class OPENGL_CALLBACK_TOGGLE_SMOOTH_SHADING:public OPENGL_CALLBACK
{
    bool smooth_shading;
    OPENGL_WORLD& world;
public:
    OPENGL_CALLBACK_TOGGLE_SMOOTH_SHADING(OPENGL_WORLD& world_input)
        :smooth_shading(false),world(world_input)
    {}

    void operator() () PHYSBAM_OVERRIDE
    {
        if(smooth_shading){
            smooth_shading=false;
            for(int i=1; i<=world.object_list.m; i++) if(world.can_toggle_smooth_shading(i))
                world.object_list(i)->Turn_Smooth_Shading_Off();}
        else{
            smooth_shading=true;
            for(int i=1; i<=world.object_list.m; i++) if(world.can_toggle_smooth_shading(i))
                world.object_list(i)->Turn_Smooth_Shading_On();}
    }
    void Print(std::ostream& out) { out << "Toggle smooth shading"; }
};
//#####################################################################
// Class OPENGL_CALLBACK_SAVE_CAMERA
//#####################################################################
class OPENGL_CALLBACK_SAVE_VIEW:public OPENGL_CALLBACK
{
    OPENGL_WORLD &world;
    std::string filename;
    bool verbose;
public:
    OPENGL_CALLBACK_SAVE_VIEW(OPENGL_WORLD &world_input, const std::string& filename_input,bool verbose_input=false)
        :world(world_input),filename(filename_input),verbose(verbose_input)
    {}

    void operator() () PHYSBAM_OVERRIDE
    {
        world.Save_View(filename, verbose);
    }
    void Print(std::ostream& out) { out << "Save view"; }
};
//#####################################################################
// Class OPENGL_CALLBACK_LOAD_VIEW
//#####################################################################
class OPENGL_CALLBACK_LOAD_VIEW:public OPENGL_CALLBACK
{
    OPENGL_WORLD &world;
    std::string filename;
    bool verbose;
public:
    OPENGL_CALLBACK_LOAD_VIEW(OPENGL_WORLD &world_input,const std::string& filename_input,bool verbose_input=false)
        :world(world_input),filename(filename_input),verbose(verbose_input)
    {}

    void operator() () PHYSBAM_OVERRIDE
    {
        world.Load_View(filename, verbose);
    }
    void Print(std::ostream& out) { out << "Load view"; }
};
//#####################################################################
// Class OPENGL_CALLBACK_STANDARD_SIZE
//#####################################################################
// GL-specific
class OPENGL_CALLBACK_STANDARD_SIZE:public OPENGL_CALLBACK
{
    OPENGL_WORLD &world;
public:

    OPENGL_CALLBACK_STANDARD_SIZE(OPENGL_WORLD &world)
        :world(world)
    {}

    void operator() () PHYSBAM_OVERRIDE
    {
        world.Resize_Window(640,480);
    }
    void Print(std::ostream& out) { out << "Resize Window to 640x480"; }
};
//#####################################################################
}
#endif

