//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GENERIC_PARSER
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/GENERIC_PARSER.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <fstream>
using namespace PhysBAM;
//#####################################################################
// Function Reset_Parser
//#####################################################################
template<class T> void GENERIC_PARSER<T>::
Reset_Parser()
{
    line=lines.begin();partial_line=*line;
}
//#####################################################################
// Function Get_Statement
//#####################################################################
template<class T> bool GENERIC_PARSER<T>::
Get_Statement(std::string& identifier,PARAMETER_LIST& list)
{
    identifier="";std::string rest;std::string data;
    for(;;){
        std::string::size_type start=partial_line.find("{");std::string::size_type end;
        if(start==std::string::npos){
            identifier+=partial_line;++line;
            if(line==lines.end()){STRING_UTILITIES::Strip_Whitespace(identifier);if(identifier=="") return false;else goto parse_error;}
            partial_line=*line;continue;}
        else {identifier+=partial_line.substr(0,start);rest=partial_line.substr(start+1);}
        //accumulate data until we find closing brace to end block
        while((end=rest.find("}"))==std::string::npos){data+=rest+"\n";++line;if(line==lines.end())goto parse_error;else rest=*line;}
        data+=rest.substr(0,end);partial_line=rest.substr(end+1);
        //LOG::cout<<"Got ident="<<identifier<<" data="<<data<<std::endl;
        std::istringstream string_stream(data.c_str());
        list.Read(string_stream);//,std::istringstream:in);
        STRING_UTILITIES::Strip_Whitespace(identifier);
        return true;}
  parse_error:
    LOG::cerr<<"Parse Error..."<<std::endl;
    return false;
}
//#####################################################################
// Function Preprocess_File
//#####################################################################
template<class T> void GENERIC_PARSER<T>::
Preprocess_File(std::string raw_filename,const int frame)
{
    // strip comments and do include files
    std::string filename=FILE_UTILITIES::Get_Frame_Filename(raw_filename,frame);
    std::ifstream input_stream(filename.c_str());
    std::string directory=FILE_UTILITIES::Get_Base_Directory_Name(filename);
    if(!input_stream) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Could not open scene file %s",filename.c_str()));
    std::string line,statement;
    while(std::getline(input_stream,line)){
        std::string::size_type comment=line.find("//");
        statement=(comment==std::string::npos)?line:line.substr(0,comment);
        STRING_UTILITIES::Strip_Whitespace(statement);
        if(statement.substr(0,9)=="#include "){
            std::string::size_type start=statement.find("\""),end=statement.rfind("\"");
            if(start==std::string::npos||end==std::string::npos)
                LOG::cerr<<filename<<" line "<<line<<": Could not open included file "<<filename<<std::endl;
            else Preprocess_File(directory+"/"+statement.substr(start+1,end-start-1),frame);}
        else lines.Append(statement);}
    input_stream.close();
}
//#####################################################################
template class GENERIC_PARSER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GENERIC_PARSER<double>;
#endif
#endif
