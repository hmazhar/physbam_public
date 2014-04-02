//#####################################################################
// Copyright 2003-2004, Eran Guendelman, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARAMETER_LIST
//##################################################################### 
// This can be used to handle parameters stored in text files and supplied in commandline.
// Usage (for combination with commandline arguments):
//
//      PARAMETER_LIST parameter_list;
//      parameter_list.Begin_Parse(argc,argv);
//      VECTOR<double,3> direction=parameter_list.Get_Parameter_And_Set_If_Undefined("light_direction,ld",VECTOR<double,3>(1,0,0),"light direction");
//      ...
//      parameter_list.End_Parse();
//      // Now parameter_list.Get_Extra_Parameter(i) will give you the i'th extra parameter appearing in the commandline
//
// Note that the compound name "light_direction,ld" in this case means that the long name of the parameter is 'light_direction' and the short name
// is 'ld'.  Either short or long names can be used in commandline arguments.  Only long name can be used in the parameter file.
// (i.e. in the commandline you can use "-ld 0 1 0" or "-light_direction 0 1 0", etc.
// End_Parse() is where some usage/error messages may be printed if something went wrong with the parsing.
//
// Commandline usage:
//      <program_name> [extraargs...] [-param1 <param1_value>] [-param2 <param2_value>] [-- extraargs...]
//
// Parameters are identified by being of the form -<name> with <name> being a non-number.
// Basically in order for commandline parsing to be less ambiguous, extra arguments can come either before the first parameter
// or after a special delimiter '--' at the end of the commandline.
//##################################################################### 
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __PARAMETER_LIST__
#define __PARAMETER_LIST__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Log/LOG.h>
namespace PhysBAM{

class PARAMETER_LIST:public NONCOPYABLE
{
private:
    class PARAMETER_INFO
    {
    public:
        std::string value_string,description;
        std::string default_value_string,short_name; // Only really useful when combined with commandline parsing
        PARAMETER_INFO() {}
        PARAMETER_INFO(const std::string& value_string,const std::string& description)
            : value_string(value_string),description(description) {}
        PARAMETER_INFO(const std::string& value_string,const std::string& description,const std::string& default_value_string,const std::string& short_name)
            : value_string(value_string),description(description),default_value_string(default_value_string),short_name(short_name) {}
    };

    std::string program_name;
    std::string usage;

    std::string first_error_parameter;
    ARRAY<PAIR<std::string,PARAMETER_INFO> > parameter_list;
    bool verbose;

    // To support combining with command line arguments
    HASHTABLE<std::string,std::string> commandline_parameter_map;
    ARRAY<std::string> extra_commandline_parameters;

public:
    PARAMETER_LIST()
        :verbose(false)
    {}

    //
    // Basic functions for adding parameters and getting their values
    //
    template<class T> void Set_Parameter(const std::string& name,const T& value,const std::string& description="")
    {std::string value_string;Value_To_String(value,value_string);
    Insert_Parameter(name,PARAMETER_INFO(value_string,description),true);}

    template<class T> T Get_Parameter(const std::string& compound_name,const T& default_value=T(),const std::string& description="",bool set_if_undefined=false)
    {std::string long_name,short_name;Get_Long_And_Short_Names(compound_name,long_name,short_name);
    std::string default_value_string;Value_To_String(default_value,default_value_string);

    // Transfer value from commandline parameters if it's there (overwriting one already in table)
    std::string value_string;
    if(Find_And_Remove_From_Commandline_Parameter_Map(long_name,short_name,value_string)){
        if(Find_Parameter(long_name)) LOG::cout<<"Note: '"<<long_name<<"' overridden by commandline value ("<<value_string<<")"<<std::endl;
        // In some cases (e.g. bools) we allow different syntax in the commandline, so we convert it here
        T value;if(String_To_Value<T>(value_string,value,true)) Value_To_String<T>(value,value_string);
        Insert_Parameter(long_name,PARAMETER_INFO(value_string,description),true);}

    // Now deal with the main parameter list as usual
    T value;int i=Find_Parameter(long_name);
    if(i){
        // Initialize stored default value and short name
        parameter_list(i).y.default_value_string=default_value_string;
        parameter_list(i).y.short_name=short_name;
        if(!String_To_Value<T>(parameter_list(i).y.value_string,value)){
            if(first_error_parameter.empty()) first_error_parameter=long_name;
            return default_value;}
        else{
            if(verbose) LOG::cout<<"[param] Got "<<long_name<<"="<<parameter_list(i).y.value_string<<std::endl;
            return value;}}
    else{
        if(set_if_undefined) Insert_Parameter(long_name,PARAMETER_INFO(default_value_string,description,default_value_string,short_name),true);
        return default_value;}
    }

    // treat char* as std::string
    std::string Get_Parameter(const std::string& compound_name,const char* default_value,const std::string& description="",bool set_if_undefined=false)
    {return Get_Parameter(compound_name,std::string(default_value),description,set_if_undefined);}

    template<class T> void Get_Parameter_In_Place(const std::string& compound_name,T& parameter,const std::string& description="",bool set_if_undefined=false)
    {parameter=Get_Parameter(compound_name,parameter,description,set_if_undefined);}

    template<class T> T Get_Parameter_And_Set_If_Undefined(const std::string& compound_name,const T& default_value=T(),const std::string& description="")
    {return Get_Parameter(compound_name,default_value,description,true);}

    int Number_Of_Extra_Parameters()
    {return extra_commandline_parameters.m;}

    template<class T> T Get_Extra_Parameter(int index)
    {if(index>extra_commandline_parameters.m){Print_Usage();exit(1);}
    T value;String_To_Value<T>(extra_commandline_parameters(index),value);return value;}

    void Set_Usage(const std::string& usage_input)
    {usage=usage_input;}

    void Set_Verbose(const bool verbose_input=true)
    {verbose=verbose_input;}

//#####################################################################
private:
    //
    // Parsing functions (define these for each type you want to support)
    //
    template<class T> static void Value_To_String(const T& value,std::string& value_string);
    template<class T> static bool String_To_Value(const std::string& value_string,T& value,bool commandline_style=false);   // returns false if parsing errors

public:
    void Begin_Parse();
    void Begin_Parse(const std::string& parameters_filename);
    void Begin_Parse(int argc,char *argv[],const std::string& default_extra_parameters_filename="");
    bool End_Parse(bool exit_on_error=true);

public:
    //
    // Utilities
    //
    void Clear();

    bool Is_Defined(std::string name) const;    // parameter exists
    void Rename(const std::string& old_name,const std::string& new_name);

    // stream should be a regular (ascii) stream
    void Write(std::ostream& output_stream) const;
    bool Read(std::istream& input_stream,bool overwrite_existing_values=false);

    // convenience functions
    void Write(const std::string& filename) const;
    bool Read(const std::string& filename,bool overwrite_existing_values=false);

private:
    void Print_Usage() const;
    void Parse_Commandline_Arguments(int argc,char *argv[]);
    int Find_Parameter(const std::string& name) const;
    void Insert_Parameter(const std::string& name,const PARAMETER_INFO& parameter_info,bool overwrite_existing_values=false);
    static void Get_Long_And_Short_Names(const std::string& compound_name,std::string& long_name,std::string& short_name);
    bool Find_And_Remove_From_Commandline_Parameter_Map(const std::string& long_name,const std::string& short_name,std::string& value_string);
};
}
#endif
#endif
