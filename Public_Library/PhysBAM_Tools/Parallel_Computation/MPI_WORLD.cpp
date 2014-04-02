//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_WORLD
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>

#ifdef USE_MPI

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <cstdlib>
#include <stdexcept>
#include <mpi.h>
#include <stdio.h>
#ifdef LAM_MPI
#include <lam_config.h>
#endif
namespace PhysBAM{
//#####################################################################
// Function MPI_Error_Handler
//#####################################################################
static int global_rank=0;
static void MPI_Error_Handler(MPI::Comm& comm,int* return_code,...) PHYSBAM_UNUSED;
static void MPI_Error_Handler(MPI::Comm& comm,int* return_code,...)
{
    MPI::Exception exception(*return_code);
    LOG::cerr<<"******************************** MPI ERROR ********************************"<<std::endl;
    LOG::cerr<<exception.Get_error_string()<<", global rank "<<global_rank<<std::endl;
    PROCESS_UTILITIES::Backtrace();
    throw std::runtime_error("MPI Error");
}
//#####################################################################
// Constructor
//#####################################################################
MPI_WORLD::
MPI_WORLD()
{
    Initialize(false);
}
//#####################################################################
// Constructor
//#####################################################################
MPI_WORLD::
MPI_WORLD(int& argc,char**& argv)
{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    Initialize(PARSE_ARGS::Find_And_Remove("-mpi",argc,argv));
#else
    Initialize(false);
#endif
}
//#####################################################################
// Function Initialize
//#####################################################################
bool Find_Version(FILE* info,const char* format,char* version)
{
    char buffer[256];
    while(fgets(buffer,sizeof(buffer)-1,info))
        if(sscanf(buffer,format,version))
            return true;
    return false;
}
void MPI_WORLD::
Initialize(bool force_mpi)
{
    bool use_mpi=force_mpi;
    // Since popen tends to muck up gdb, we provide a way to skip the mpi version check
    bool check_version=!getenv("PHYSBAM_SKIP_MPI_VERSION_CHECK");
#ifdef LAM_MPI
    if(getenv("LAMRANK")){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"Detected LAM Runtime Environment"<<std::endl;
#endif
        use_mpi=true;
        if(check_version){
            FILE* laminfo=popen("laminfo -version lam full","r");
            if(!laminfo) PHYSBAM_FATAL_ERROR("LAM/MPI version check failed: couldn't run laminfo");
            char version[21];
            int length=fscanf(laminfo," LAM/MPI: %20s",version);
            if(!length) PHYSBAM_FATAL_ERROR("LAM/MPI version check failed: couldn't parse laminfo output");
            if(strncmp(version,LAM_VERSION,20)) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("LAM/MPI version check failed: compiled with %s, run with %s",LAM_VERSION,version));
            pclose(laminfo);}}
#elif defined(OPEN_MPI)
    if(getenv("LAMRANK"))
        PHYSBAM_FATAL_ERROR("MPI version check failed: compiled with OpenMPI, run in LAM/MPI");
    if(getenv("OMPI_MCA_universe") || getenv("OMPI_UNIVERSE_SIZE")){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"Detected OpenMPI Runtime Environment"<<std::endl;
#endif
        use_mpi=true;
        if(check_version){
            std::string compiled_version=STRING_UTILITIES::string_sprintf("%d.%d.%d",OMPI_MAJOR_VERSION,OMPI_MINOR_VERSION,OMPI_RELEASE_VERSION);
            FILE* ompi_info=popen("ompi_info -version ompi full","r");
            if(!ompi_info) PHYSBAM_FATAL_ERROR("Open MPI version check failed: couldn't run ompi_info");
            char version[21];
            if(!Find_Version(ompi_info," Open MPI: %20s",version))
                PHYSBAM_FATAL_ERROR("Open MPI version check failed: couldn't parse ompi_info output");
            if(compiled_version!=version) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Open MPI version check failed: compiled with %s, run with %s",compiled_version.c_str(),version));
            pclose(ompi_info);}}
#else
//#error Unrecognized MPI package
#endif
    if(!use_mpi){initialized=false;return;}
    int argc=0;char** argv=0; // these are ignored by lam and openmpi
    int status=MPI_Init(&argc,&argv);
    if(status!=MPI_SUCCESS) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("MPI_Init failed: status %d",status));
    initialized=true;
    global_rank=rank=MPI::COMM_WORLD.Get_rank();

#if defined(OPEN_MPI) && OMPI_MAJOR_VERSION>1 || (OMPI_MAJOR_VERSION==1 && OMPI_MINOR_VERSION>=3)
    // This works with recent version of OpenMPI, but used to break with earlier versions of LAM/OpenMPI (not sure which)
    MPI::Errhandler handler=MPI::COMM_WORLD.Create_errhandler(MPI_Error_Handler);
    MPI::COMM_WORLD.Set_errhandler(handler);
    handler.Free();
#endif
}
//#####################################################################
// Destructor
//#####################################################################
MPI_WORLD::
~MPI_WORLD()
{
    if(initialized) MPI_Finalize();
}
//#####################################################################
// Function Initialized
//#####################################################################
bool MPI_WORLD::
Initialized()
{
    return MPI::Is_initialized();
}
//#####################################################################
}
#else
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
MPI_WORLD::
MPI_WORLD(int& argc,char**& argv)
    :initialized(false),rank(0)
{
    LOG::cerr<<"Not compiled with USE_MPI."<<std::endl;
}
//#####################################################################
// Destructor
//#####################################################################
MPI_WORLD::
~MPI_WORLD()
{}
//#####################################################################
// Function Initialized
//#####################################################################
bool MPI_WORLD::
Initialized()
{
    return false;
}
//#####################################################################
}
#endif
