//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// EXAMPLE
//#####################################################################
template<class TV> EXAMPLE<TV>::
EXAMPLE(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),last_frame(120),frame_rate(24),frame_title(""),write_substeps_level(-1),write_first_frame(true),write_last_frame(true),write_time(true),
    output_directory("output"),data_directory("../../Public_Data/Archives"),auto_restart(false),restart(false),restart_frame(0),write_output_files(true),write_frame_title(true),
    abort_when_dt_below(0),parse_args(0),mpi_world(0),want_mpi_world(false),need_finish_logging(false),test_number(0),fixed_dt((T)0),substeps_delay_frame(-1),substeps_delay_level(-1)
{
    ARRAY<std::string> directory_tokens;
    STRING_UTILITIES::Split(FILE_UTILITIES::Get_Working_Directory(),"/",directory_tokens);
    bool is_archived_project=false;
    for(int i=1;i<=directory_tokens.Size();++i) if(STRING_UTILITIES::Ends_With(directory_tokens(i), "Archives")) is_archived_project=true;
    if(is_archived_project) data_directory = "../" + data_directory;
}
//#####################################################################
// ~EXAMPLE
//#####################################################################
template<class TV> EXAMPLE<TV>::
~EXAMPLE()
{
    delete mpi_world;
    if(need_finish_logging) LOG::Finish_Logging();
    delete parse_args;
}
template<class TV> typename TV::SCALAR EXAMPLE<TV>::
Time_At_Frame(const int frame) const
{
    return initial_time+(frame-first_frame)/frame_rate;
}
template<class TV> void EXAMPLE<TV>::
Clamp_Time_Step_With_Target_Time(const T time,const T target_time,T& dt,bool& done,const T min_dt,bool* min_dt_failed)
{
    if(dt<min_dt){dt=min_dt;if(min_dt_failed) *min_dt_failed=true;}
    if(time+dt>=target_time){dt=target_time-time;done=true;}
    else if(time+2*dt>=target_time) dt=min(dt,(T).51*(target_time-time));
}
template<class TV> void EXAMPLE<TV>::
Set_Write_Substeps_Level(const int level)
{
    write_substeps_level=level;
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(level);
}
template<class TV> void EXAMPLE<TV>::
Write_Frame_Title(const int frame) const
{
    if(write_frame_title) FILE_UTILITIES::Write_To_Text_File(STRING_UTILITIES::string_sprintf("%s/%d/frame_title",output_directory.c_str(),frame),frame_title);
}
template<class TV> void EXAMPLE<TV>::
Limit_Dt(T& dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Log_Parameters() const
{
    LOG::SCOPE scope("EXAMPLE parameters");
    LOG::cout<<"initial_time="<<initial_time<<std::endl;
    LOG::cout<<"first_frame="<<first_frame<<std::endl;
    LOG::cout<<"last_frame="<<last_frame<<std::endl;
    LOG::cout<<"frame_rate="<<frame_rate<<std::endl;
    LOG::cout<<"auto_restart="<<auto_restart<<std::endl;
    LOG::cout<<"restart="<<restart<<std::endl;
    LOG::cout<<"restart_frame="<<restart_frame<<std::endl;
    LOG::cout<<"write_output_files="<<write_output_files<<std::endl;
    LOG::cout<<"write_first_frame="<<write_first_frame<<std::endl;
    LOG::cout<<"write_last_frame="<<write_last_frame<<std::endl;
    LOG::cout<<"write_time="<<write_time<<std::endl;
    LOG::cout<<"write_frame_title="<<write_frame_title<<std::endl;
    LOG::cout<<"write_substeps_level="<<write_substeps_level<<std::endl;
    LOG::cout<<"output_directory="<<output_directory<<std::endl;
    LOG::cout<<"data_directory="<<data_directory<<std::endl;
    LOG::cout<<"frame_title="<<frame_title<<std::endl;
    LOG::cout<<"abort_when_dt_below="<<abort_when_dt_below<<std::endl;
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Register_Options()
{
    if(!parse_args) return;
    parse_args->Set_Extra_Arguments(-1,"<example number>");
    parse_args->Add_String_Argument("-o","","output directory");
    parse_args->Add_String_Argument("-d","","data directory");
    parse_args->Add_Integer_Argument("-restart",0,"frame","restart frame");
    parse_args->Add_Option_Argument("-auto_restart","restart from last_frame");
    parse_args->Add_Integer_Argument("-substeps",-1,"level","substep output level");
    parse_args->Add_Integer_Argument("-delay_substeps",-1,"frame","delay substeps until later frame");
    parse_args->Add_Integer_Argument("-first_frame",0,"frame","first frame");
    parse_args->Add_Integer_Argument("-last_frame",0,"frame","last frame");
    parse_args->Add_Double_Argument("-framerate",24,"frame rate");
    parse_args->Add_Option_Argument("-query_output","print the output directory and exit");
    parse_args->Add_Integer_Argument("-v",1<<30,"level","verbosity level");
    parse_args->Add_Option_Argument("-nolog","disable log.txt");
    parse_args->Add_Double_Argument("-dt",0,"fix the time step size to this value.");
    if(mpi_world) parse_args->Add_Option_Argument("-all_verbose","all mpi processes write to stdout (not just the first)");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Parse_Options()
{
    if(!parse_args) return;
    test_number=Subexample(1);
    if(parse_args->Is_Value_Set("-d")) data_directory=parse_args->Get_String_Value("-d");
    else if(const char* d=getenv("PHYSBAM_DATA_DIRECTORY")) data_directory=d;

    int verbosity=parse_args->Get_Integer_Value("-v");
    if(mpi_world && !parse_args->Is_Value_Set("-all_verbose") && mpi_world->initialized && mpi_world->rank) verbosity=0;
    need_finish_logging=true;
    LOG::Initialize_Logging(verbosity<10,false,verbosity,!parse_args->Is_Value_Set("-nolog"));
}
//#####################################################################
// Function Override_Options
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Override_Options()
{
    if(parse_args->Is_Value_Set("-o")) output_directory=parse_args->Get_String_Value("-o");
    if(parse_args->Is_Value_Set("-restart")){restart=true;restart_frame=parse_args->Get_Integer_Value("-restart");}
    if(parse_args->Is_Value_Set("-auto_restart")) auto_restart=true;
    if(parse_args->Is_Value_Set("-delay_substeps")){
        substeps_delay_frame=parse_args->Get_Integer_Value("-delay_substeps");
        substeps_delay_level=parse_args->Get_Integer_Value("-substeps");}
    else if(parse_args->Is_Value_Set("-substeps")) Set_Write_Substeps_Level(parse_args->Get_Integer_Value("-substeps"));
    if(parse_args->Is_Value_Set("-first_frame")) first_frame=parse_args->Get_Integer_Value("-first_frame");
    if(parse_args->Is_Value_Set("-framerate")) frame_rate=parse_args->Get_Double_Value("-framerate");
    if(parse_args->Is_Value_Set("-query_output")){LOG::cout<<output_directory;exit(0);}
    if(parse_args->Is_Value_Set("-dt")) fixed_dt=(T)parse_args->Get_Double_Value("-dt");

    if(!parse_args->Is_Value_Set("-nolog")){
        if(!restart && !auto_restart) FILE_UTILITIES::Create_Directory(output_directory);
        FILE_UTILITIES::Create_Directory(output_directory+"/common");
        LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",restart);}
}
//#####################################################################
// Function Parse_Late_Options
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Parse_Late_Options()
{
    if(!parse_args) return;
    if(parse_args->Is_Value_Set("-last_frame")) last_frame=parse_args->Get_Integer_Value("-last_frame");
}
//#####################################################################
// Function Parse
//#####################################################################
template<class TV> void EXAMPLE<TV>::
Parse(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);

    if(want_mpi_world) mpi_world=new MPI_WORLD(argc,argv);
    parse_args=new PARSE_ARGS;
    Register_Options();

    parse_args->Parse(argc,argv);
    std::string print_args=parse_args->Print_Arguments(argc,argv);

    Parse_Options();
    LOG::cout<<print_args<<std::endl;
    Override_Options();
}
//#####################################################################
// Function Subexample
//#####################################################################
template<class TV> int EXAMPLE<TV>::
Subexample(const int default_example) const
{
    if(parse_args->Num_Extra_Args()<1) return default_example;
    int parsed_value;
    if(STRING_UTILITIES::String_To_Value(parse_args->Extra_Arg(1),parsed_value)) return parsed_value;
    throw VALUE_ERROR("The argument is not an integer.");
}
//#####################################################################
template class EXAMPLE<VECTOR<float,1> >;
template class EXAMPLE<VECTOR<float,2> >;
template class EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EXAMPLE<VECTOR<double,1> >;
template class EXAMPLE<VECTOR<double,2> >;
template class EXAMPLE<VECTOR<double,3> >;
#endif
