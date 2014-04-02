//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __DRIVER__
#define __DRIVER__
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{


template<class TV> class EXAMPLE;

template<class TV>
class DRIVER
{
    typedef typename TV::SCALAR T;
protected:
    T time;
public:
    EXAMPLE<TV>& example;
    int current_frame;
    int output_number;

    DRIVER(EXAMPLE<TV>& example);
    virtual ~DRIVER();

    void Write_First_Frame(const int frame) const
    {if(example.write_first_frame && frame==example.first_frame) FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");}

    void Write_Last_Frame(const int frame) const
    {if(example.write_last_frame) FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");}

    void Write_Time(const int frame) const
    {if(example.write_time) FILE_UTILITIES::Write_To_File(example.stream_type,STRING_UTILITIES::string_sprintf("%s/%d/time",example.output_directory.c_str(),frame),time);}

    T Time() const
    {return time;}

    void Set_Time(const T time_input) 
    {time=time_input;}
    
//#####################################################################
    virtual void Execute_Main_Program();
    virtual void Initialize();
    virtual void Advance_To_Target_Time(const T target_time)=0;
    virtual void Write_Substep(const std::string& title,const int substep,const int level=0);
    virtual void Read_Time(const int frame);
    virtual void Read_Last_Frame();
    virtual void Simulate_To_Frame(const int frame);
    virtual void Write_Output_Files(const int frame);
//#####################################################################
};
}
#endif
