//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/DRIVER.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void*writer,const std::string& title,int substep,int level)
{
    ((DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// DRIVER
//#####################################################################
template<class TV> DRIVER<TV>::
DRIVER(EXAMPLE<TV>& example)
    :example(example),output_number(example.first_frame)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// ~DRIVER
//#####################################################################
template<class TV> DRIVER<TV>::
~DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void DRIVER<TV>::
Initialize()
{
    // setup time
    if(example.auto_restart){Read_Last_Frame();example.restart=true;}
    if(example.restart){current_frame=example.restart_frame;Read_Time(current_frame);}
    else current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        LOG::cout<<"Writing substep ["<<example.frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;
        Write_Output_Files(++output_number);example.frame_title="";}
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(current_frame<frame){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Write_Output_Files(++output_number);
        current_frame++;}
}
//#####################################################################
// Function Read_Time
//#####################################################################
template<class TV> void DRIVER<TV>::
Read_Time(const int frame)
{
    time=example.Time_At_Frame(frame);
    std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/time",example.output_directory.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){
        T corrected_time;
        FILE_UTILITIES::Read_From_File(example.stream_type,filename,corrected_time);
        if(abs(time-corrected_time)>(T)1e-4*abs(time)){ // only adjust time if significantly different from default in order to get deterministic restarts
            time=corrected_time;
            // adjust initial time so that Simulate_To_Frame() returns correct time (essential when writing substeps)
            example.initial_time=time-(frame-example.first_frame)/example.frame_rate;}}
}
//#####################################################################
// Function Read_Last_Frame
//#####################################################################
template<class TV> void DRIVER<TV>::
Read_Last_Frame()
{
    std::string filename=STRING_UTILITIES::string_sprintf("%s/common/last_frame",example.output_directory.c_str());
    if(FILE_UTILITIES::File_Exists(filename))
        FILE_UTILITIES::Read_From_Text_File(filename, example.restart_frame);
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    Write_First_Frame(frame);
    example.Write_Output_Files(frame);
    Write_Time(frame);
    Write_Last_Frame(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),example.frame_title);
}
//#####################################################################
template class DRIVER<VECTOR<float,1> >;
template class DRIVER<VECTOR<float,2> >;
template class DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DRIVER<VECTOR<double,1> >;
template class DRIVER<VECTOR<double,2> >;
template class DRIVER<VECTOR<double,3> >;
#endif
