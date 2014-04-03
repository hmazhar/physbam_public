//#####################################################################
// Class OPENGL_WORLD
//##################################################################### 
// The philosophy: provide an easy way to visualize scenes (that may change
// under the user's control, but are not interactive per se --- this is
// primarily for data visualization).
// Create an OPENGL_WORLD object, add objects and lights as you like,
// bind keys to change the default (e.g. to allow for switching between
// files), and then call Run_Visualization().
// You can also call Initialize() and Start_Main_Loop() instead of
// Run_Visualization() if you need to do some things after the OpenGL
// context has been created and before the main loop is started.
//
// Prompt mechanism:
//   Call Prompt_User with the prompt message and a callback.
//   User will type in response, and either press ENTER or ESC.  
//   Callback will then be called.
//   If prompt was finished with ESC (i.e. it was aborted) then
//   prompt_response will be null.  If it was finished with ENTER
//   then prompt_response will be the response string.
//#####################################################################
#ifndef __OPENGL_WORLD__
#define __OPENGL_WORLD__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_KEY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <ctime>
namespace PhysBAM{

class OPENGL_WINDOW;
class OPENGL_CALLBACK;
class OPENGL_SELECTION;
class OPENGL_LIGHT;
class OPENGL_MOUSE_HANDLER;
template<class T> class OPENGL_ARCBALL;

class OPENGL_KEY_BINDING_CATEGORY
{
public:
    std::string name;
    int priority;
    ARRAY<PAIR<OPENGL_KEY,OPENGL_CALLBACK*> > key_bindings;

    OPENGL_KEY_BINDING_CATEGORY() {}
    OPENGL_KEY_BINDING_CATEGORY(const std::string &name,int priority) : name(name),priority(priority) {}
};

class OPENGL_WORLD
{
    typedef float T;
    typedef VECTOR<T,3> TV;
public:
    enum FILL_MODE { DRAW_FILLED, DRAW_WIREFRAME, DRAW_FILLED_AND_WIREFRAME };

    typedef void (*PROCESS_HITS_CB)(GLint hits, GLuint buffer[]);

private:
    bool initialized;
public:

    // objects
    ARRAY<OPENGL_OBJECT*> object_list; // does not own objects
    ARRAY<bool> can_toggle_smooth_shading;
    ARRAY<bool> use_bounding_box;

    // lights
    OPENGL_COLOR ambient_light;
    ARRAY<OPENGL_LIGHT*> lights;

    // window parameters
    float fovy;
    bool mode_2d;
    bool load_names_for_selection;
    OPENGL_WINDOW* window;
private:

    // display parameters
    bool stereo_mode;
    float stereo_offset;
    int fill_mode;
    bool enable_lighting_for_wireframe;
    bool white_background;
    ARRAY<PLANE<float>*> clipping_planes;

    // string displays
    ARRAY<std::string> strings_to_print;
    bool display_strings;
    bool show_object_names;
    bool display_object_names_in_corner; // In display name mode, display them in the corner
    bool view_auto_help;

    // Idle/Timer Data
    OPENGL_CALLBACK* idle_callback;
    int timer_id;
    float idle_delay;
    float idle_timer;
    float view_target_timer;
    float frame_counter_timer;
    int frames_rendered,frames_per_second;
    bool show_frames_per_second;

    // Camera
    bool left_handed_coordinate_system;
    float nearclip_factor,farclip_factor,nearclip,farclip;
    OPENGL_ARCBALL<float>* arcball; // This maintains the arcballs info
    float camera_distance;
    TV target_position;
    MATRIX<float,4> arcball_matrix; // This is used to extract the last extracted matrix from the arcball
    MATRIX<float,4> rotation_matrix; // This stores the current rotation of the ball

    // Mouse and Camera Interaction
    int zoom_direction,translation_direction;
    int oldmousex,oldmousey;
    bool do_mouse_rotation;
    bool do_mouse_zoom;
    bool do_mouse_target_xy, do_mouse_target_z;
    TV target_x_drag_vector;
    TV target_y_drag_vector;
    OPENGL_MOUSE_HANDLER* external_mouse_handler;
    bool shift_was_pressed; 
    bool ctrl_was_pressed; 

    // Keyboard Interaction
    ARRAY<ARRAY<OPENGL_CALLBACK*>,VECTOR<int,2> > key_bindings;
    std::string current_key_binding_category;
    int current_key_binding_category_priority;
    ARRAY<OPENGL_KEY_BINDING_CATEGORY> key_bindings_by_category;

    // Prompting
    bool prompt_mode;
    std::string prompt;
public:
    std::string prompt_response;
    bool prompt_response_success;
private:
    OPENGL_CALLBACK* prompt_response_cb;

    // Selection
    PROCESS_HITS_CB process_hits_cb;
    bool selection_mode;  // selection stuff
    TV* current_selection; // pointer to current selection item    

public:
    OPENGL_WORLD();
    ~OPENGL_WORLD();
    static OPENGL_WORLD* Singleton();
    void Run_Visualization(const std::string& window_title="OpenGL Visualization");
    void Initialize(const std::string& window_title="OpenGL Visualization",const int width=640,const int height=480,const bool offscreen=false);
    void Initialize_Glut_Independent();
    void Clear_All_Objects();
    void Add_Object(OPENGL_OBJECT* object,bool include_bounding_box=true,bool toggle_smooth_shading=false);
    void Clear_All_Lights();
    void Set_Ambient_Light(float value);
    void Set_Ambient_Light(const OPENGL_COLOR& color);
    void Add_Light(OPENGL_LIGHT* light);
    void Set_Key_Binding_Category(const std::string &category);
    void Set_Key_Binding_Category_Priority(int priority=1);

    void Bind_Key(const OPENGL_KEY& key,OPENGL_CALLBACK* callback);
    void Bind_Key(const std::string& key,OPENGL_CALLBACK* callback);
    void Append_Bind_Key(const OPENGL_KEY& key,OPENGL_CALLBACK* callback);
    void Append_Bind_Key(const std::string& key,OPENGL_CALLBACK* callback);
    void Unbind_Key(const OPENGL_KEY& key);
    void Unbind_Keys(const std::string& keys);

    void Set_Idle_Callback(OPENGL_CALLBACK* callback,const float delay);
    void Set_View_Target_Timer(const float view_target_timer_input);
private:
    void Prepare_For_Idle();
    void Handle_Idle();
    void Handle_Timer();
public:

    void Clear_Strings();
    void Add_String(const std::string& s);
    void Set_Zoom_Direction(bool up_is_zoom_in);
    void Set_Translation_Direction(bool up_is_move_in);
    void Set_Lighting_For_Wireframe(bool enable_flag);
    void Set_2D_Mode(bool mode=true);
    void Set_Process_Hits_Callback(PROCESS_HITS_CB process_hits_cb_input);

    TV Get_Camera_Position();
    TV Get_Target_Position();
    void Get_View_Frame(TV &view_forward, TV &view_up, TV &view_right);
    void Set_View_Frame(const TV &view_forward, const TV &view_up, const TV &view_right);
    void Get_Look_At(TV &camera, TV &target, TV &up);
    void Set_Look_At(const TV &camera, const TV &target, const TV &up);
    void Save_View(const std::string& filename, bool verbose = false);
    bool Load_View(const std::string& filename, bool verbose = false);
    RAY<VECTOR<float,3> > Ray_Through_Normalized_Image_Coordinate(VECTOR<float,2> coordinates);
    void Set_Left_Handed(const bool left_handed=false);

    RANGE<TV> Scene_Bounding_Box();
    void Center_Camera_On_Bounding_Box(const RANGE<TV>& bounding_box,const bool adjust_distance);
    void Reset_Camera_Orientation(const bool reset_up_vector_only=false);
    void Center_Camera_On_Scene();

    VECTOR<float,2> Convert_Mouse_Coordinates(int x, int y);

    void Handle_Display_Prompt_Only();
    void Render_World(bool selecting,bool swap_buffers=true);
    void Handle_Reshape_Main();
    void Handle_Keypress_Main(const OPENGL_KEY& key,int x,int y);
    void Handle_Keypress_Prompt(unsigned char raw_key);
    void Handle_Click_Main(int button,int state,int x,int y,bool ctrl_pressed,bool shift_pressed);
    void Handle_Drag_Main(int x,int y);
    void Display_Target(const int in_color=1);
    void Display_Auto_Help();
    static void Draw_Transparent_Text_Box(const ARRAY<std::string> &strings,const VECTOR<int,2> &top_left_corner,int vspace,void* font,const OPENGL_COLOR &color);
    void Display_Strings(const ARRAY<std::string> &strings,const OPENGL_COLOR &color=OPENGL_COLOR::White(),bool draw_transparent_box=true,int horizontal_offset=0,int vspace=18,void* font=GLUT_BITMAP_HELVETICA_18);
    void Display_Strings(bool draw_transparent_box=true);
    void Display_Object_Names();
    void Display_Object_Names_In_Corner();
    static int Get_Number_Of_Valid_Hits(GLint hits,GLuint buffer[],int buffer_size);
    static void Print_Hits(GLint hits, GLuint buffer[]);
    void Prepare_For_Target_XY_Drag();

    void Set_External_Mouse_Handler(OPENGL_MOUSE_HANDLER* mouse_handler=0);

    void Update_Clipping_Planes();
    GLenum Add_Clipping_Plane(const PLANE<float> &plane);
    void Set_Clipping_Plane(GLenum id,const PLANE<float> &plane);
    void Remove_Clipping_Plane(GLenum id);
    void Remove_All_Clipping_Planes();

    // To use this you need to set load_names_for_selection
    void Get_Selections(ARRAY<OPENGL_SELECTION*> &selections, GLint hits, GLuint buffer[]);

    void Save_Screen(const std::string& filename,const bool use_back_buffer,int jpeg_quality=90);
    template<int d> void Get_Image(ARRAY<VECTOR<T,d>,VECTOR<int,2> >& image,const bool use_back_buffer);
    static void Resize_Window(const int width,const int height);

    void Display_Prompt_Strings();
    void Prompt_User(const std::string& prompt,OPENGL_CALLBACK* prompt_response_cb){Prompt_User(prompt,prompt_response_cb,"");} // workaround for MSVC's inability to handle defaults
    void Prompt_User(const std::string& prompt,OPENGL_CALLBACK* prompt_response_cb,const std::string& default_response);

    // Basic callbacks
    void Toggle_Background();
    void Toggle_Help();
    void Toggle_Show_Frames_Per_Second();
    void Resize_To_Standard_Size();
    void Quit();
    DEFINE_CALLBACK_CREATOR(OPENGL_WORLD,Toggle_Background);
    DEFINE_CALLBACK_CREATOR(OPENGL_WORLD,Toggle_Help);
    DEFINE_CALLBACK_CREATOR(OPENGL_WORLD,Toggle_Show_Frames_Per_Second);
    DEFINE_CALLBACK_CREATOR(OPENGL_WORLD,Resize_To_Standard_Size);
    DEFINE_CALLBACK_CREATOR(OPENGL_WORLD,Quit);
public:
    friend class OPENGL_WINDOW;
    friend class OPENGL_WINDOW_GLUT;
//#####################################################################
};
}
#endif

