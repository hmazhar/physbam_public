//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Images/MOV_FILE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <climits>
#include <cstdlib>
#include <fstream>
using namespace PhysBAM;

ANIMATED_VISUALIZATION::ANIMATED_VISUALIZATION()
    :animation_enabled(true),play(false),loop(false),fixed_frame_rate(false),start_frame(0),stop_frame(INT_MAX),frame(0),frame_rate(24),frame_increment(1),last_frame_filename("")
{
    if(MOV_WRITER<float>::Enabled()) saved_frame_filename="capture.mov";
    else if(IMAGE<float>::Is_Supported(".png")) saved_frame_filename="capture.%05d.png";
    else if(IMAGE<float>::Is_Supported(".jpg")) saved_frame_filename="capture.%05d.jpg";
    else saved_frame_filename="capture.%05d.ppm";
}

void ANIMATED_VISUALIZATION::
Add_Arguments(PARSE_ARGS& parse_args)
{
    BASIC_VISUALIZATION::Add_Arguments(parse_args);

    parse_args.Add_Integer_Argument("-jpeg_quality",95,"jpeg quality settings");
    parse_args.Add_String_Argument("-so",saved_frame_filename,"save frames output");
    parse_args.Add_Integer_Argument("-start_frame",0,"start frame");
    parse_args.Add_Integer_Argument("-stop_frame",INT_MAX,"stop frame");
    parse_args.Add_Integer_Argument("-fps",24,"frames per second");
}

void ANIMATED_VISUALIZATION::
Parse_Arguments(PARSE_ARGS& parse_args)
{
    BASIC_VISUALIZATION::Parse_Arguments(parse_args);

    jpeg_quality=parse_args.Get_Integer_Value("-jpeg_quality");
    if(parse_args.Is_Value_Set("-so")) saved_frame_filename=parse_args.Get_String_Value("-so");
    if(parse_args.Is_Value_Set("-start_frame")) start_frame=parse_args.Get_Integer_Value("-start_frame");
    if(parse_args.Is_Value_Set("-stop_frame")) stop_frame=parse_args.Get_Integer_Value("-stop_frame");
    if(parse_args.Is_Value_Set("-fps")) frame_rate=parse_args.Get_Integer_Value("-fps");
}

void ANIMATED_VISUALIZATION::
Add_OpenGL_Initialization()
{
    BASIC_VISUALIZATION::Add_OpenGL_Initialization();

    Update_OpenGL_Strings();
}

void ANIMATED_VISUALIZATION::
Initialize_Components_And_Key_Bindings()
{
    BASIC_VISUALIZATION::Initialize_Components_And_Key_Bindings();

    opengl_world.Set_Key_Binding_Category("Default Keys (ANIMATED_VISUALIZATION)");
    opengl_world.Set_Key_Binding_Category_Priority(100);

    opengl_world.Bind_Key('p',Toggle_Play_CB("Play/Pause"));
    opengl_world.Bind_Key('P',Toggle_Loop_CB("Loop"));
    opengl_world.Bind_Key('r',Reset_CB("Reset"));
    opengl_world.Bind_Key('s',Next_Frame_CB("Next frame"));
    opengl_world.Bind_Key('S',Next_Frame_CB("Next frame"));
    opengl_world.Bind_Key("^s",Prev_Frame_CB("Prev frame"));
    opengl_world.Bind_Key('g',Goto_Frame_CB("Goto frame"));
    opengl_world.Bind_Key('f',Toggle_Fixed_Frame_Rate_CB("Toggle fixed frame rate"));
    opengl_world.Bind_Key('z',Goto_Last_Frame_CB("Goto last frame"));
    opengl_world.Bind_Key("^d",Capture_Frames_CB("Capture frames (to images)"));

    opengl_world.Set_Key_Binding_Category("User-Defined Keys");
    opengl_world.Set_Key_Binding_Category_Priority(1);
}

bool ANIMATED_VISUALIZATION::
Valid_Frame(int frame_input)
{
    if(last_frame_filename!=""){
        std::ifstream last_frame_file(last_frame_filename.c_str());
        if(last_frame_file){
            int last_frame;
            if((last_frame_file>>last_frame) && frame_input>last_frame) return false;}}
    // It's valid if it's valid for any of our animated components...
    for(int i=1;i<=component_list.m;i++)
        if(component_list(i)->Is_Animated() && component_list(i)->Valid_Frame(frame_input)) return true;
    return false;
}

void ANIMATED_VISUALIZATION::
Goto_Start_Frame()
{
    Set_Frame(start_frame);
}

void ANIMATED_VISUALIZATION::
Goto_Last_Frame()
{
    if(last_frame_filename!=""){
        int last_frame;
        std::ifstream last_frame_file(last_frame_filename.c_str());
        if(last_frame_file && (last_frame_file>>last_frame))
            Set_Frame(last_frame);}
}

void ANIMATED_VISUALIZATION::
Render_Offscreen()
{
    Capture_Frames(saved_frame_filename,start_frame,stop_frame,jpeg_quality,false);
}

OPENGL_EPS_OUTPUT<float>* PhysBAM::opengl_eps_output = 0;

void ANIMATED_VISUALIZATION::
Capture_Frames(const std::string& filename_pattern,int capture_start_frame,int capture_end_frame,int jpeg_quality,bool swap_buffers)
{
    bool use_eps=STRING_UTILITIES::IEnds_With(filename_pattern,".eps");
    bool movie=STRING_UTILITIES::IEnds_With(filename_pattern,".mov");
    MOV_WRITER<float>* mov=0;
    if(movie && MOV_WRITER<float>::Enabled()){
        mov=new MOV_WRITER<float>(filename_pattern,24);
        LOG::cout<<"Capturing to quicktime '"<<filename_pattern<<"'"<<std::endl;}

    LOG::cout<<"Capturing frames "<<capture_start_frame<<" to ";
    if(capture_end_frame==INT_MAX) LOG::cout<<"last valid frame";else LOG::cout<<capture_end_frame;
    LOG::cout<<" into '"<<filename_pattern<<"'"<<std::endl;
    Set_Frame(capture_start_frame);
    for(;;){
        if(use_eps){
            std::string filename=STRING_UTILITIES::string_sprintf(filename_pattern.c_str(),frame);
            opengl_eps_output = new OPENGL_EPS_OUTPUT<float>(filename);}
        opengl_world.Render_World(false,swap_buffers);
        glFinish();
        if(mov){
            LOG::cout<<"  Frame "<<frame<<std::endl;
            ARRAY<VECTOR<float,3> ,VECTOR<int,2> > image;
            opengl_world.Get_Image(image,swap_buffers);
            if(mov) mov->Add_Frame(image);}
        else if(!use_eps){
            std::string filename=STRING_UTILITIES::string_sprintf(filename_pattern.c_str(),frame);
            LOG::cout<<"Capturing frame "<<frame<<" to "<<filename<<std::endl;
            opengl_world.Save_Screen(filename,swap_buffers,jpeg_quality);}
        else{delete opengl_eps_output;opengl_eps_output=0;}
        if(!animation_enabled || frame==capture_end_frame || !Valid_Frame(frame+1)) break;
        Set_Frame(frame+1);}
    delete mov;
}

//#####################################################################
// CALLBACKS
//#####################################################################

void ANIMATED_VISUALIZATION::
Set_Frame(int frame_input)
{
    if(!animation_enabled) return;
    frame=frame_input;

    Pre_Frame_Extra();

    for(int i=1;i<=component_list.m;i++){
#ifdef NDEBUG
        int attempts=0;
        bool done=false;
        while(!done && attempts++ < 10){
            try{
                component_list(i)->Set_Frame(frame);
                done=true;}
            catch(std::exception& error){
                LOG::cerr<<"Read error: "<<error.what()<<std::endl;
                LOG::cerr<<"Retrying in 1/10 s..."<<std::endl;
                PROCESS_UTILITIES::Sleep(.1);}}
#else
        // don't do try-catch in debug to allow gdb to catch it
        component_list(i)->Set_Frame(frame);
#endif
    }

    // Update/Invalidate selections
    for(int i=1;i<=component_list.m;i++){
        bool delete_selection=false;
        OPENGL_SELECTION* new_selection=component_list(i)->Create_Or_Destroy_Selection_After_Frame_Change(current_selection,delete_selection);
        if(new_selection){Set_Current_Selection(new_selection);break;}
        else if(delete_selection){Set_Current_Selection(0);break;}}

    Set_Frame_Extra();
    Update_OpenGL_Strings();
}

void ANIMATED_VISUALIZATION::
Update_OpenGL_Strings()
{
    opengl_world.Clear_Strings();
    if(animation_enabled) opengl_world.Add_String(STRING_UTILITIES::string_sprintf("frame %d",frame)+(frame_title.empty()?"":": "+frame_title));
    BASIC_VISUALIZATION::Update_OpenGL_Strings();
}

void ANIMATED_VISUALIZATION::
Next_Frame()
{
    if(!animation_enabled) return;
    bool valid=Valid_Frame(frame+frame_increment);
    if(valid || loop){
        Set_Frame(valid?frame+frame_increment:start_frame);
        opengl_world.Set_Idle_Callback(play?Next_Frame_CB():0,fixed_frame_rate?(float)1/frame_rate:0);}
    else if(play) opengl_world.Set_Idle_Callback(Next_Frame_CB(),.2);
}

void ANIMATED_VISUALIZATION::
Prev_Frame()
{
    if(!animation_enabled) return;
    if(Valid_Frame(frame-frame_increment))
        Set_Frame(frame-frame_increment);
}

void ANIMATED_VISUALIZATION::
Goto_Frame_Prompt()
{
    if(!opengl_world.prompt_response.empty()){
        int input_frame;
        STRING_UTILITIES::String_To_Value(opengl_world.prompt_response,input_frame);
        if(Valid_Frame(input_frame)) Set_Frame(input_frame);}
}

void ANIMATED_VISUALIZATION::
Goto_Frame()
{
    if(!animation_enabled) return;
    opengl_world.Prompt_User("Goto frame: ",Goto_Frame_Prompt_CB(),"");
}

void ANIMATED_VISUALIZATION::
Reset()
{
    if(!animation_enabled) return;
    Set_Frame(start_frame);
    if(play) Toggle_Play(); // Stop playing
}

void ANIMATED_VISUALIZATION::
Toggle_Play()
{
    if(!animation_enabled) return;
    play=!play;
    opengl_world.Set_Idle_Callback(play?Next_Frame_CB():0,fixed_frame_rate?(float)1/frame_rate:0);
}

void ANIMATED_VISUALIZATION::
Toggle_Loop()
{
    loop=!loop;
}

void ANIMATED_VISUALIZATION::
Toggle_Fixed_Frame_Rate()
{
    fixed_frame_rate=!fixed_frame_rate;
    if(fixed_frame_rate) last_frame_time=0;
    opengl_world.Set_Idle_Callback(play?Next_Frame_CB():0,fixed_frame_rate?(float)1/frame_rate:0);
}

void ANIMATED_VISUALIZATION::
Capture_Frames_Prompt()
{
    if(!opengl_world.prompt_response_success) return;

    bool done=false;
    if(capture_frames_prompt_state.step==1){
        if(!opengl_world.prompt_response.empty()) capture_frames_prompt_state.filename_pattern=opengl_world.prompt_response;
        capture_frames_prompt_state.step=2;
        opengl_world.Prompt_User(STRING_UTILITIES::string_sprintf("Start frame [%d]: ",start_frame),Capture_Frames_Prompt_CB(),"");}
    else if(capture_frames_prompt_state.step==2){
        if(!opengl_world.prompt_response.empty()) STRING_UTILITIES::String_To_Value(opengl_world.prompt_response,capture_frames_prompt_state.start_frame);
        capture_frames_prompt_state.step=3;
        opengl_world.Prompt_User("End frame [last valid frame]: ",Capture_Frames_Prompt_CB(),"");}
    else if(capture_frames_prompt_state.step==3){
        if(!opengl_world.prompt_response.empty()) STRING_UTILITIES::String_To_Value(opengl_world.prompt_response,capture_frames_prompt_state.end_frame);
        capture_frames_prompt_state.step=4;
        if(STRING_UTILITIES::toupper(FILE_UTILITIES::Get_File_Extension(capture_frames_prompt_state.filename_pattern))=="JPG")
            opengl_world.Prompt_User(STRING_UTILITIES::string_sprintf("JPEG quality [%d]: ",jpeg_quality),Capture_Frames_Prompt_CB(),"");
        else done=true;}
    else if(capture_frames_prompt_state.step==4){
        if(opengl_world.prompt_response.empty()) STRING_UTILITIES::String_To_Value(opengl_world.prompt_response,capture_frames_prompt_state.jpeg_quality);
        done=true;}

    if(done){
        Capture_Frames(capture_frames_prompt_state.filename_pattern,capture_frames_prompt_state.start_frame,capture_frames_prompt_state.end_frame,
            capture_frames_prompt_state.jpeg_quality);
        if(system(STRING_UTILITIES::string_sprintf("pbp %s &",capture_frames_prompt_state.filename_pattern.c_str()).c_str())){};}
}

void ANIMATED_VISUALIZATION::
Capture_Frames()
{
    // Fill with defaults
    capture_frames_prompt_state.filename_pattern=saved_frame_filename;
    capture_frames_prompt_state.start_frame=start_frame;
    capture_frames_prompt_state.end_frame=INT_MAX;
    capture_frames_prompt_state.jpeg_quality=jpeg_quality;
    capture_frames_prompt_state.step=1;
    opengl_world.Prompt_User("Capture filename [" + saved_frame_filename + "]: ",Capture_Frames_Prompt_CB(),"");
}
