//#####################################################################
// Class OPENGL_UNICALLBACK
//##################################################################### 
#ifndef __OPENGL_UNICALLBACK__
#define __OPENGL_UNICALLBACK__

#ifdef USE_FLTK
#include <PhysBAM_Tools/Log/LOG.h>

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
    class classname##func : public PhysBAM::OPENGL_UNICALLBACK { \
        public: \
            classname##func(classname *obj, const std::string &help_string) : obj(obj), help_string(help_string) {} \
            void operator() () PHYSBAM_OVERRIDE { obj->func(); } \
            void Print(std::ostream &out) { out << help_string; } \
        private: \
             classname *obj; \
             std::string help_string; }; \
    PhysBAM::OPENGL_UNICALLBACK *func##_CB(const std::string &help_string = "Unknown") \
    { return new classname##func(this, help_string); }

namespace PhysBAM{

class OPENGL_UNICALLBACK
{
public:
    virtual ~OPENGL_UNICALLBACK() {}
    virtual void operator() () = 0;
    virtual void Print(std::ostream& out) { out << "Unknown"; }
};

//#####################################################################
// now a few useful callbacks that are used in the default bindings

//#####################################################################
// Class OPENGL_UNICALLBACK_QUIT
//#####################################################################
class OPENGL_UNICALLBACK_QUIT:public OPENGL_UNICALLBACK
{
public:
    void operator() () PHYSBAM_OVERRIDE {exit(0);}
    void Print(std::ostream& out) { out << "Quit"; }
};
//#####################################################################
// Class OPENGL_UNICALLBACK_TOGGLE
//#####################################################################
class OPENGL_UNICALLBACK_TOGGLE:public OPENGL_UNICALLBACK
{
    bool *state;
    std::string help_string;
public:
    OPENGL_UNICALLBACK_TOGGLE(bool *state_input, const std::string &help_string_input="Toggle state") 
        : state(state_input), help_string(help_string_input) {}
    void operator() () PHYSBAM_OVERRIDE {*state=!*state;}
    void Print(std::ostream& out) { out << help_string << " (" << (*state?"1":"0") << ")"; }
};
//#####################################################################
// Class OPENGL_UNICALLBACK_CYCLE
//#####################################################################
class OPENGL_UNICALLBACK_CYCLE:public OPENGL_UNICALLBACK
{
    int *state;
    int min_value, max_value;
    std::string help_string;
public:
    OPENGL_UNICALLBACK_CYCLE(int *state, int min_value, int max_value, const std::string &help_string_input="Cycle state") 
        : state(state), min_value(min_value), max_value(max_value), help_string(help_string_input) {}
    void operator() () PHYSBAM_OVERRIDE { if (++(*state) > max_value) *state = min_value; }
    void Print(std::ostream& out) { out << help_string << " (" << *state << ")"; }
};
//#####################################################################
// Class OPENGL_UNICALLBACK_SET_VALUE
//#####################################################################
template<class T>
class OPENGL_UNICALLBACK_SET_VALUE:public OPENGL_UNICALLBACK
{
    T *state, value;
    std::string help_string;
public:
    OPENGL_UNICALLBACK_SET_VALUE(T *state_input, T value_input, const std::string &help_string_input="Set value") 
        : state(state_input),value(value_input),help_string(help_string_input) {}
    void operator() () PHYSBAM_OVERRIDE {*state=value;}
    void Print(std::ostream& out) { out << help_string << " (" << *state << "=" << value << ")"; }
};
//#####################################################################
// Class OPENGL_UNICALLBACK_ZOOM
//#####################################################################
// actually a misnomer - this is panning in/out, not zooming in the true sense of the word
class OPENGL_UNICALLBACK_ZOOM:public OPENGL_UNICALLBACK
{
    const float factor;
    float *camera_distance, *nearclip, *farclip;
public:
    OPENGL_UNICALLBACK_ZOOM(const float factor_input,float *camera_distance_input,float *nearclip_input,float *farclip_input)
        : factor(factor_input),camera_distance(camera_distance_input),nearclip(nearclip_input),farclip(farclip_input)
    {}

    void operator() () PHYSBAM_OVERRIDE {*camera_distance*=factor;*nearclip*=factor;*farclip*=factor;}
    void Print(std::ostream& out) { out << "Zoom by " << factor; }
};
//#####################################################################
// Class OPENGL_UNICALLBACK_WHITE_BACKGROUND
//#####################################################################
// GL-specific
class OPENGL_UNICALLBACK_WHITE_BACKGROUND:public OPENGL_UNICALLBACK
{
    bool white;
public:
    OPENGL_UNICALLBACK_WHITE_BACKGROUND()
        :white(false)
    {}

    void operator() () PHYSBAM_OVERRIDE
    {
        white=!white;
        if(white){
            glClearColor(1.0,1.0,1.0,1.0);
//          OPENGL_WORLD::wireframe_color=OPENGL_COLOR::Black();
        }
        else{
            glClearColor(0.0,0.0,0.0,0.0);
//          OPENGL_WORLD::wireframe_color=OPENGL_COLOR::Gray(0.7);
        }
    }
    void Print(std::ostream& out) { out << "Toggle white background"; }
};

//=====================================================================
// OPENGL_WORLD specific callbacks
//=====================================================================

//#####################################################################
// Class OPENGL_UNICALLBACK_SAVE_SCREEN
//#####################################################################
class OPENGL_UNICALLBACK_SAVE_SCREEN:public OPENGL_UNICALLBACK
{
    OPENGL_SUBSPACE &world;
    const std::string basename;
    int counter;
public:
    OPENGL_UNICALLBACK_SAVE_SCREEN(OPENGL_SUBSPACE &world_input,const std::string& basename_input="opengl",int counter_start=0)
        : world(world_input),basename(basename_input),counter(counter_start)
    {}

    void operator() () PHYSBAM_OVERRIDE
    {
        std::string defaultname,filename;
        defaultname=STRING_UTILITIES::string_sprintf("%s%03d.ppm",basename.c_str(),counter);
        LOG::cout<<"Filename for saving image in .ppm or .jpg (if supported) formats (return for default, - for cancel)\n";
        LOG::cout<<"[%s]: "<<defaultname;
        getline(std::cin,filename);
        size_t len=filename.length();if(filename[len-1]=='\n'){--len;filename=filename.substr(0,len);}
        if(len==0) filename=defaultname;
        else if(filename=="-") return;
        else{} // should parse user input to make a better default basename next time
        LOG::cout<<"saving to %s..."<<filename;
//      world.Save_Screen(filename);
        LOG::cout<<"\n";
        ++counter;
    }
    void Print(std::ostream& out) { out << "Save screen"; }
};
//#####################################################################
// Class OPENGL_UNICALLBACK_MOVE_TARGET
//#####################################################################
class OPENGL_UNICALLBACK_MOVE_TARGET:public OPENGL_UNICALLBACK
{
    OPENGL_SUBSPACE &world;
    const VECTOR<float,3> motion_step;
    const float *motion_factor;
    VECTOR<float,3> *target;
    float *view_target_timer;
public:
    OPENGL_UNICALLBACK_MOVE_TARGET(OPENGL_SUBSPACE &world_input, const VECTOR<float,3>& motion_step_input,const float *motion_factor_input,
                                VECTOR<float,3> *target_input, float *view_target_timer_input)
        : world(world_input),motion_step(motion_step_input),motion_factor(motion_factor_input),target(target_input),
          view_target_timer(view_target_timer_input)
    {}

    void operator() ()  PHYSBAM_OVERRIDE
    {
        VECTOR<float,3> view_forward, view_up, view_right;
        world.Get_View_Frame(view_forward, view_up, view_right);
        (*target) += *motion_factor * (motion_step.x * view_right + motion_step.y * view_up + motion_step.z * view_forward);
        *view_target_timer=1;
    }
    void Print(std::ostream& out) { out << "Move target " << motion_step; }
};
//#####################################################################
// Class OPENGL_UNICALLBACK_TOGGLE_SMOOTH_SHADING
//#####################################################################
class OPENGL_UNICALLBACK_TOGGLE_SMOOTH_SHADING:public OPENGL_UNICALLBACK
{
    bool smooth_shading;
    OPENGL_SUBSPACE& world;
public:
    OPENGL_UNICALLBACK_TOGGLE_SMOOTH_SHADING(OPENGL_SUBSPACE& world_input)
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
// Class OPENGL_UNICALLBACK_SAVE_CAMERA
//#####################################################################
class OPENGL_UNICALLBACK_SAVE_VIEW:public OPENGL_UNICALLBACK
{
    OPENGL_SUBSPACE &world;
    char filename[512];
    bool verbose;
public:
    OPENGL_UNICALLBACK_SAVE_VIEW(OPENGL_SUBSPACE &world_input, const char *filename_input, bool verbose_input=false)
        : world(world_input), verbose(verbose_input)
    { strcpy(filename, filename_input); }

    void operator() () PHYSBAM_OVERRIDE
    {
        world.Save_View(filename, verbose);
    }
    void Print(std::ostream& out) { out << "Save view"; }
};
//#####################################################################
// Class OPENGL_UNICALLBACK_LOAD_VIEW
//#####################################################################
class OPENGL_UNICALLBACK_LOAD_VIEW:public OPENGL_UNICALLBACK
{
    OPENGL_SUBSPACE &world;
    char filename[512];
    bool verbose;
public:
    OPENGL_UNICALLBACK_LOAD_VIEW(OPENGL_SUBSPACE &world_input, const char *filename_input, bool verbose_input=false)
        : world(world_input), verbose(verbose_input)
    { strcpy(filename, filename_input); }

    void operator() () PHYSBAM_OVERRIDE
    {
        world.Load_View(filename, verbose);
    }
    void Print(std::ostream& out) { out << "Load view"; }
};
//#####################################################################
}
#endif
#endif
