//#####################################################################
// Copyright 2003-2004, Eran Guendelman, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARAMETER_LIST
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_ROTATION.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <cctype>
#include <fstream>
#include <sstream>
namespace PhysBAM{
//#####################################################################
// Utility functions (in unnamed namespace)
//#####################################################################
namespace
{
    template<class T> static void streamable_value_to_string(const T& value,std::string& value_string)
    {std::ostringstream value_string_stream;
    value_string_stream<<value;
    value_string=value_string_stream.str();}

    template<class T> static bool streamable_string_to_value(const std::string& value_string,T& value)
    {std::istringstream value_string_stream(value_string);
    return value_string_stream>>value;}

    static bool is_dash_parameter_name(const std::string &str,std::string *parameter_name=0)
    {if(str.length()>1 && str[0]=='-'){
        std::string tmp_parameter_name=STRING_UTILITIES::Stripped_Whitespace(str.substr(1));
        if(!STRING_UTILITIES::Is_Number(tmp_parameter_name) && tmp_parameter_name!="-"){
            if(parameter_name) *parameter_name=tmp_parameter_name;
            return true;}}
    return false;}

    static bool is_end_of_parameters_delimiter(const std::string &str)
    {return str=="--";}
}

//#####################################################################
// Parsing for streamable types
//   Basically these are types which already have << and >> defined
//#####################################################################
#define DEFINE_STREAMABLE_PARAMETER(T) \
    template<> void PARAMETER_LIST::Value_To_String(const T& value,std::string& value_string) \
    {streamable_value_to_string(value,value_string);} \
    template<> bool PARAMETER_LIST::String_To_Value(const std::string& value_string,T& value,bool commandline_style) \
    {return streamable_string_to_value(value_string,value);}

#define DEFINE_STREAMABLE_VECTOR_PARAMETER(T,d) \
    template<> void PARAMETER_LIST::Value_To_String(const VECTOR<T,d>& value,std::string& value_string) \
    {streamable_value_to_string(value,value_string);} \
    template<> bool PARAMETER_LIST::String_To_Value(const std::string& value_string,VECTOR<T,d>& value,bool commandline_style) \
    {return streamable_string_to_value(value_string,value);}
#define COMMA ,
DEFINE_STREAMABLE_PARAMETER(double);
DEFINE_STREAMABLE_PARAMETER(float);
DEFINE_STREAMABLE_PARAMETER(int);
DEFINE_STREAMABLE_PARAMETER(QUATERNION<float>);
DEFINE_STREAMABLE_PARAMETER(QUATERNION<double>);
typedef VECTOR<float,3> FV3;DEFINE_STREAMABLE_PARAMETER(ROTATION<FV3>);
typedef VECTOR<double,3> FD3;DEFINE_STREAMABLE_PARAMETER(ROTATION<FD3>);
DEFINE_STREAMABLE_PARAMETER(RANGE<VECTOR<int COMMA 3> >);
DEFINE_STREAMABLE_PARAMETER(RANGE<VECTOR<float COMMA 3> >);
DEFINE_STREAMABLE_PARAMETER(RANGE<VECTOR<double COMMA 3> >);
DEFINE_STREAMABLE_VECTOR_PARAMETER(double,2);
DEFINE_STREAMABLE_VECTOR_PARAMETER(double,3);
DEFINE_STREAMABLE_VECTOR_PARAMETER(float,2);
DEFINE_STREAMABLE_VECTOR_PARAMETER(float,3);
DEFINE_STREAMABLE_VECTOR_PARAMETER(int,2);
DEFINE_STREAMABLE_VECTOR_PARAMETER(int,3);
//#####################################################################
// Parsing for bools
//   Both "true"/"false" and "1"/"0" are supported
//   On the commandline, the existence of the parameter is enough to imply "true"
//#####################################################################
template<> void PARAMETER_LIST::
Value_To_String(const bool &value,std::string &value_string)
{
    value_string=(value) ? "true" : "false";
}
template<> bool PARAMETER_LIST::
String_To_Value(const std::string &value_string,bool &value,bool commandline_style)
{
    std::string value_string_stripped=STRING_UTILITIES::Stripped_Whitespace(value_string);
    if((commandline_style && value_string_stripped=="") || value_string_stripped=="true" || value_string_stripped=="1") value=true;
    else if(value_string_stripped=="false" || value_string_stripped=="0") value=false;
    else return false;
    return true;
}
//#####################################################################
// Parsing for strings
//#####################################################################
template<> void PARAMETER_LIST::
Value_To_String(const std::string& value,std::string& value_string)
{
    value_string="\"";value_string.append(value);value_string.append("\"");
}
template<> bool PARAMETER_LIST::
String_To_Value(const std::string& value_string,std::string& value,bool commandline_style)
{
    if(value_string[0]=='"'){
        std::string::size_type second_quotes=value_string.find('"',1);
        if(second_quotes==std::string::npos){ LOG::cerr<<"Could not find matching quotes in string"<<std::endl;return false;}
        value=value_string.substr(1,second_quotes-1);
        return true;}
    else{
        std::istringstream value_string_stream(value_string);
        value_string_stream>>value;
        return true;}
}
//#####################################################################
// Parsing for other types
//   Create specialization of PARAMETER_LIST_PARSER and instantiate using DEFINE_PARSER_PARAMETER
//#####################################################################
template<class T>
struct PARAMETER_LIST_PARSER
{
    static void Value_To_String(const T& value,std::string& value_string);
    static bool String_To_Value(const std::string& value_string,T& value,bool commandline_style);
};
#define DEFINE_PARSER_PARAMETER(...) \
    template<> void PARAMETER_LIST::Value_To_String(const __VA_ARGS__& value,std::string& value_string) \
    {PARAMETER_LIST_PARSER<__VA_ARGS__>::Value_To_String(value,value_string);} \
    template<> bool PARAMETER_LIST::String_To_Value(const std::string& value_string,__VA_ARGS__& value,bool commandline_style) \
    {return PARAMETER_LIST_PARSER<__VA_ARGS__>::String_To_Value(value_string,value,commandline_style);}
//#####################################################################
// Parsing for RANGE<VECTOR<T,2> >
//#####################################################################
template<class T> struct PARAMETER_LIST_PARSER<RANGE<VECTOR<T,2> > >
{
    static void Value_To_String(const RANGE<VECTOR<T,2> >& value,std::string& value_string)
    {
        std::ostringstream value_string_stream;value_string_stream<<value.min_corner.x<<" "<<value.max_corner.x<<" "<<value.min_corner.y<<" "<<value.max_corner.y;
        value_string=value_string_stream.str();
    }
    static bool String_To_Value(const std::string& value_string,RANGE<VECTOR<T,2> >& value,bool commandline_style)
    {
        std::istringstream value_string_stream(value_string);value_string_stream>>value.min_corner.x>>value.max_corner.x>>value.min_corner.y>>value.max_corner.y;
        return value_string_stream!=0;
    }
};
DEFINE_PARSER_PARAMETER(RANGE<VECTOR<int,2> >);
DEFINE_PARSER_PARAMETER(RANGE<VECTOR<float,2> >);
DEFINE_PARSER_PARAMETER(RANGE<VECTOR<double,2> >);
//#####################################################################
// Function Begin_Parse
//#####################################################################
void PARAMETER_LIST::
Begin_Parse()
{
    Clear();
}
//#####################################################################
// Function Begin_Parse
//#####################################################################
void PARAMETER_LIST::
Begin_Parse(const std::string &parameters_filename)
{
    Clear();
    Read(parameters_filename);
}
//#####################################################################
// Function Begin_Parse
//#####################################################################
void PARAMETER_LIST::
Begin_Parse(int argc,char *argv[],const std::string &default_extra_parameters_filename)
{
    Clear();
    Parse_Commandline_Arguments(argc,argv);

    // Read in file specified by -param
    bool param_defined=Is_Defined("param");
    std::string param_filename=Get_Parameter_And_Set_If_Undefined<std::string>("param",default_extra_parameters_filename,"read extra parameters from");
    if(!param_filename.empty() && (param_defined || FILE_UTILITIES::File_Exists(param_filename))){
        LOG::cout<<"Note: Reading parameters from "<<param_filename<<std::endl;
        Read(param_filename);}
}
//#####################################################################
// Function End_Parse
//#####################################################################
bool PARAMETER_LIST::
End_Parse(bool exit_on_error)
{
    // Print help if requested
    if(Get_Parameter<bool>("help",false,"print help")){Print_Usage();exit(1);}

    // Check for parsing errors
    bool has_error=(!first_error_parameter.empty() || commandline_parameter_map.Size());
    if(has_error) Print_Usage();
    if(!first_error_parameter.empty()){
        LOG::cout<<"ERROR: Problems parsing parameter '"<<first_error_parameter<<"'"<<std::endl;
        int i=Find_Parameter(first_error_parameter);
        if(i) LOG::cout<<"       (Got value string '"<<parameter_list(i).y.value_string<<"')"<<std::endl;}
    if(commandline_parameter_map.Size())
        LOG::cout<<"ERROR: Unrecognized commandline parameters:";
    LOG::cout<<commandline_parameter_map<<std::endl;
    if(has_error && exit_on_error) exit(1);
    return has_error;
}
//#####################################################################
// Function Clear
//#####################################################################
void PARAMETER_LIST::
Clear()
{
    parameter_list.Remove_All();
    commandline_parameter_map.Clean_Memory();
    extra_commandline_parameters.Remove_All();
    first_error_parameter="";
}
//#####################################################################
// Function Is_Defined
//#####################################################################
bool PARAMETER_LIST::
Is_Defined(std::string name) const
{
    // Either defined in commandline arguments or parameter list itself
    return commandline_parameter_map.Contains(name) || Find_Parameter(name);
}
//#####################################################################
// Function Rename
//#####################################################################
void PARAMETER_LIST::
Rename(const std::string &old_name,const std::string &new_name)
{
    PHYSBAM_ASSERT(Is_Defined(old_name) && !Is_Defined(new_name));
    int i=Find_Parameter(old_name);
    if(i) parameter_list(i).x=new_name;
}
//#####################################################################
// Function Write (stream)
//#####################################################################
void PARAMETER_LIST::
Write(std::ostream& output_stream) const
{
    for(int i=1;i<=parameter_list.m;i++){
        output_stream<<parameter_list(i).x<<" = "<<parameter_list(i).y.value_string;
        if(parameter_list(i).y.description.length()>0)
            output_stream<<" // "<<parameter_list(i).y.description;
        output_stream<<std::endl;}
}
//#####################################################################
// Function Read (stream)
//#####################################################################
bool PARAMETER_LIST::
Read(std::istream& input_stream,bool overwrite_existing_values)
{
    std::string line;
    while(std::getline(input_stream,line)){
        std::string::size_type comment=line.find("//");
        std::string statement,description;
        if(comment==std::string::npos){statement=STRING_UTILITIES::Stripped_Whitespace(line);description="";}
        else{statement=STRING_UTILITIES::Stripped_Whitespace(line.substr(0,comment));description=STRING_UTILITIES::Stripped_Whitespace(line.substr(comment+2));}
        if(statement.empty())continue;
        std::string::size_type equals=statement.find("=");
        if(equals==std::string::npos){LOG::cerr<<"PARAMETER_LIST: Did not find '='"<<std::endl;return false;}
        std::string name=STRING_UTILITIES::Stripped_Whitespace(statement.substr(0,equals)),value=STRING_UTILITIES::Stripped_Whitespace(statement.substr(equals+1));
        Insert_Parameter(name,PARAMETER_INFO(value,description),overwrite_existing_values);}
    return true;
}
//#####################################################################
// Function Write (filename)
//#####################################################################
void PARAMETER_LIST::
Write(const std::string &filename) const
{
    std::ofstream output_stream(filename.c_str());
    Write(output_stream);
}
//#####################################################################
// Function Read (filename)
//#####################################################################
bool PARAMETER_LIST::
Read(const std::string &filename,bool overwrite_existing_values)
{
    std::ifstream input_stream(filename.c_str());
    if(input_stream) return Read(input_stream,overwrite_existing_values);
    else{LOG::cerr<<"Could not open parameter file "<<filename<<std::endl;return false;}
}
//#####################################################################
// Function Print_Usage
//#####################################################################
void PARAMETER_LIST::
Print_Usage() const
{
    if(!usage.empty())LOG::cout<<"Usage: "<<program_name<<" [-param parameterfile] "<<usage<<"\n";
    LOG::cout<<"Parameters:"<<std::endl;
    int width=0;
    for(int i=1;i<=parameter_list.m;i++){
        int len=(int)parameter_list(i).x.length()+1;
        if(parameter_list(i).y.short_name.length()>0) len+=(int)parameter_list(i).y.short_name.length()+4;
        if(len>width)width=len;}
    for(int i=1;i<=parameter_list.m;i++){
        std::string name;
        if(parameter_list(i).y.short_name.length()>0) name="-"+parameter_list(i).y.short_name+", ";
        name+="-"+parameter_list(i).x;
        LOG::cout.flags(std::ios::left);LOG::cout.width(width+2);
        LOG::cout<<name;
        LOG::cout<<(!parameter_list(i).y.description.empty()?parameter_list(i).y.description+" ":"");
        if(parameter_list(i).y.default_value_string.length()>0) LOG::cout<<"(default: "<<parameter_list(i).y.default_value_string<<")";
        LOG::cout<<std::endl;}
}
//#####################################################################
// Function Parse_Commandline_Arguments
//#####################################################################
void PARAMETER_LIST::
Parse_Commandline_Arguments(int argc,char *argv[])
{
    program_name=argv[0];
    std::string parameter_name;int i=1;
    // extra parameters before first named parameter
    while(i<argc && !is_dash_parameter_name(argv[i]) && !is_end_of_parameters_delimiter(argv[i])){
        extra_commandline_parameters.Append(argv[i]);i++;}
    // named parameters
    while(i<argc && is_dash_parameter_name(argv[i],&parameter_name)){
        std::string parameter_value;i++;
        while(i<argc && !is_end_of_parameters_delimiter(argv[i]) && !is_dash_parameter_name(argv[i])){
            if(!parameter_value.empty())parameter_value+=" ";
            parameter_value+=argv[i];i++;}
        commandline_parameter_map.Set(parameter_name,parameter_value);}
    // extra parameters after end_of_parameters_delimiter
    if(i<argc && is_end_of_parameters_delimiter(argv[i])){
        i++; // skip it
        while(i<argc){extra_commandline_parameters.Append(argv[i]);i++;}}
}
//#####################################################################
// Function Find_Parameter
//#####################################################################
int PARAMETER_LIST::
Find_Parameter(const std::string& name) const
{
    for(int i=1;i<=parameter_list.m;i++)if(parameter_list(i).x==name)return i;
    return 0;
}
//#####################################################################
// Function Insert_Parameter
//#####################################################################
void PARAMETER_LIST::
Insert_Parameter(const std::string& name,const PARAMETER_INFO& parameter_info,bool overwrite_existing_values)
{
    int i=Find_Parameter(name);
    if(!i) parameter_list.Append(PAIR<std::string,PARAMETER_INFO>(name,parameter_info));
    else if(overwrite_existing_values) parameter_list(i).y=parameter_info;
}
//#####################################################################
// Function Get_Long_And_Short_Names
//#####################################################################
void PARAMETER_LIST::
Get_Long_And_Short_Names(const std::string& name,std::string& long_name,std::string& short_name)
{
    std::string::size_type comma_pos=name.find(",");
    if(comma_pos!=std::string::npos){
        long_name=name.substr(0,comma_pos);
        short_name=name.substr(comma_pos+1);}
    else{
        long_name=name;
        short_name="";}
}
//#####################################################################
// Function Find_And_Remove_From_Commandline_Parameter_Map
//#####################################################################
// Look for either long or short forms in commandline parameters (long form overrides)
bool PARAMETER_LIST::
Find_And_Remove_From_Commandline_Parameter_Map(const std::string& long_name,const std::string& short_name,std::string& value_string)
{
    bool short_exists=commandline_parameter_map.Get(short_name,value_string);
    if(short_exists) commandline_parameter_map.Delete(short_name);
    bool long_exists=commandline_parameter_map.Get(long_name,value_string);
    if(long_exists) commandline_parameter_map.Delete(long_name);
    if(short_exists && long_exists)
        LOG::cerr<<"Warning: Both '"<<long_name<<"' and '"<<short_name<<"' found in commandline. ('"<<long_name<<"' will override."<<std::endl;
    return short_exists || long_exists;
}
//#####################################################################
};
#endif
