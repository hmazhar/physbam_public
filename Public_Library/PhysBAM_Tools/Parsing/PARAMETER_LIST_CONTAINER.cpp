//#####################################################################
// Copyright 2004. Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST_CONTAINER.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <fstream>
using namespace PhysBAM;
//#####################################################################
// Get_Next_Parameter_List
//#####################################################################
bool PARAMETER_LIST_CONTAINER::
Get_Next_Parameter_List(std::string& identifier,PARAMETER_LIST& list)
{
    list.Clear();
    identifier="";std::string rest;
    while(true){ // Find opening brace
        std::string::size_type start=partial_line.find("{");
        if(start==std::string::npos){
            identifier+=partial_line;++line;
            if(line==lines.end()){STRING_UTILITIES::Strip_Whitespace(identifier);if(identifier=="")goto parse_done;else goto parse_error;}
            partial_line=*line;}
        else {identifier+=partial_line.substr(0,start);rest=partial_line.substr(start+1);break;}}

    {std::string data;std::string::size_type end;
    while((end=rest.find("}"))==std::string::npos){ //accumulate data until we find closing brace to end block
        data+=rest+"\n";++line;if(line==lines.end())goto parse_error;else rest=*line;}
    data+=rest.substr(0,end);partial_line=rest.substr(end+1);
    std::istringstream string_stream(data.c_str());
    list.Read(string_stream);
    STRING_UTILITIES::Strip_Whitespace(identifier);
    return true;}

  parse_done:
    return false;

  parse_error:
    LOG::cerr<<"PARAMETER_LIST_CONTAINER: Parse Error"<<std::endl;
    return false;
}
//#####################################################################
// Pre_Process_File
//#####################################################################
// strip comments and do include files
void PARAMETER_LIST_CONTAINER::
Pre_Process_File(const std::string& filename) 
{
    std::ifstream input_stream(filename.c_str());
    std::string directory=FILE_UTILITIES::Get_Base_Directory_Name(filename);
    if(!input_stream)LOG::cerr<<"Could not open parameter list container file "<<filename<<std::endl;
    else{
        std::string line,statement;
        while(std::getline(input_stream,line)){
            std::string::size_type comment=line.find("//");
            statement=(comment==std::string::npos)?line:line.substr(0,comment);
            STRING_UTILITIES::Strip_Whitespace(statement);
            if(statement.substr(0,9)=="#include "){
                std::string::size_type start=statement.find("\""),end=statement.rfind("\"");
                if(start==std::string::npos||end==std::string::npos)
                    LOG::cerr<<filename<<" line "<<line<<": Could not open included file "<<filename<<std::endl;
                else Pre_Process_File(directory+"/"+statement.substr(start+1,end-start-1));}
            else lines.Append(statement);}
        input_stream.close();}
}
//#####################################################################
#endif
