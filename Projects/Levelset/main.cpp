#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "WATER_DRIVER.h"
#include "WATER_EXAMPLE.h"

using namespace PhysBAM;

template<class TV> void Execute_Main_Program(STREAM_TYPE& stream_type,PARSE_ARGS& parse_args,MPI_WORLD& mpi_world)
{ 
    typedef VECTOR<int,TV::dimension> TV_INT;

    WATER_EXAMPLE<TV>* example=new WATER_EXAMPLE<TV>(stream_type,parse_args.Get_Integer_Value("-threads"));

    int scale=parse_args.Get_Integer_Value("-scale");
    RANGE<TV> range(TV(),TV::All_Ones_Vector());TV_INT counts=TV_INT::All_Ones_Vector()*scale;
    example->Initialize_Grid(counts,range);
    example->restart=parse_args.Get_Integer_Value("-restart");
    example->last_frame=parse_args.Get_Integer_Value("-e");
    example->write_substeps_level=parse_args.Get_Integer_Value("-substep");
    example->write_debug_data=true;
    
    TV point1=TV::All_Ones_Vector()*(TV::dimension==2?.65:.45),point2=TV::All_Ones_Vector()*.75;point1(1)=0;point2(1)=.05;
    example->source.min_corner=point1;example->source.max_corner=point2;

    if(mpi_world.initialized){
        example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid,3);
        if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("/%d",(mpi_world.rank+1));}

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
    
    WATER_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
}

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));

    MPI_WORLD mpi_world(argc,argv);

    PARSE_ARGS parse_args;
    parse_args.Add_Integer_Argument("-restart",0,"restart frame");
    parse_args.Add_Integer_Argument("-scale",100,"fine scale grid resolution");
    parse_args.Add_Integer_Argument("-substep",-1,"output-substep level");
    parse_args.Add_Integer_Argument("-e",100,"last frame");
    parse_args.Add_Integer_Argument("-threads",1,"number of threads");
    parse_args.Add_Option_Argument("-3d","run in 3 dimensions");

    parse_args.Parse(argc,argv);
    parse_args.Print_Arguments(argc,argv);
    
    if(parse_args.Is_Value_Set("-3d")){
        Execute_Main_Program<VECTOR<T,3> >(stream_type,parse_args,mpi_world);}
    else{
        Execute_Main_Program<VECTOR<T,2> >(stream_type,parse_args,mpi_world);}

    return 0;
}
