//#####################################################################
// Copyright 2004-2008, Nipun Kwatra, Frank Losasso, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TIMER
//##################################################################### 
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <cstdio>
#include <iostream>
namespace PhysBAM{

#ifdef _WIN32 // we are running on windows
#pragma comment(lib, "winmm.lib")
#include <windows.h>
double Get_Current_Time(){__int64 time;QueryPerformanceCounter((LARGE_INTEGER*)&time);return (double)time;}
double Initialize_Timer(){__int64 frequency;QueryPerformanceFrequency((LARGE_INTEGER*)&frequency);return 1000./(double)frequency;}
#endif
#if defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__) // we are running on linux
#include <sys/resource.h>
#include <sys/time.h>
double Get_Current_Time(){timeval tv;gettimeofday(&tv,0);return tv.tv_sec+1e-6*tv.tv_usec;}
double Initialize_Timer(){return 1e3;} // to convert seconds to milliseconds
#endif

TIMER* TIMER::singleton_instance;

//#####################################################################
// Constructor
//#####################################################################
TIMER::
TIMER()
{
    timers.Resize(512);
    resolution=Initialize_Timer();
    timers(1).elapsed=timers(1).start=0;
    free_timers=IDENTITY_ARRAY<>(timers.Size());
    overhead=0;timers(1).accumulator=0;
    for(int i=1;i<=100;i++){Start(1);Stop(1);}
    overhead=timers(1).accumulator*.01;
}
//#####################################################################
// Destructor
//#####################################################################
TIMER::
~TIMER()
{}
//#####################################################################
// Function Register_Timer
//#####################################################################
int TIMER::
Register_Timer()
{
    if(!free_timers.Size()){LOG::cerr<<"No more timers available, timing information may be incorrect"<<std::endl;return 0;}
    
    int id=free_timers.Pop();

    double time=Get_Current_Time();
    timers(id).start=time;
    timers(id).elapsed=time;
    timers(id).accumulator=0;
    return id;
}
//#####################################################################
// Function Release_Timer
//#####################################################################
void TIMER::
Release_Timer(const int id)
{
    free_timers.Append(id);
}
//#####################################################################
// Function Get_Time
//#####################################################################
double TIMER::
Get_Time()
{
    return Get_Current_Time()*resolution;
}
//#####################################################################
// Function Get_Total_Time_Since_Registration
//#####################################################################
double TIMER::
Get_Total_Time_Since_Registration(const int id)
{
    return (Get_Current_Time()-timers(id).start)*resolution;
}
//#####################################################################
// Function Peek_And_Reset_Time
//#####################################################################
double TIMER::
Peek_And_Reset_Time(const int id)
{
    double time=Get_Current_Time();
    double time_elapsed=(time-timers(id).elapsed)*resolution;
    timers(id).elapsed=time;
    return time_elapsed;
}
//#####################################################################
// Function Reset_Time
//#####################################################################
void TIMER::
Reset_Time(const int id)
{
    double time=Get_Current_Time();
    timers(id).elapsed=time;
}
//#####################################################################
// Function Peek_Time
//#####################################################################
double TIMER::
Peek_Time(const int id)
{
    return (Get_Current_Time()-timers(id).elapsed)*resolution;
}
//#####################################################################
// Function Start
//#####################################################################
void TIMER::
Start(const int id)
{
    timers(id).start=Get_Current_Time();
}
//#####################################################################
// Function Stop
//#####################################################################
void TIMER::
Stop(const int id)
{
    double time=Get_Current_Time();
    timers(id).accumulator+=time-timers(id).start-overhead;
}
//#####################################################################
// Function Print_Stats
//#####################################################################
void TIMER::
Print_Stats(const int id,const char* str)
{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    LOG::cout<<"Total elapsed time for "<<str<<" is "<<timers(id).accumulator*resolution*.001<<" s."<<std::endl;
#endif
}
//#####################################################################
}
