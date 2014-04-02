//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace STRING_UTILITIES
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <cctype>
#include <cstdio>
namespace PhysBAM{
namespace STRING_UTILITIES{

std::string Stripped_Whitespace(const std::string& str)
{int first=0,last=(int)str.length()-1;
while(first<=last && isspace(str[first]))first++;
while(last>=0 && isspace(str[last]))last--;
return str.substr(first,last-first+1);}

bool Is_Number(const std::string& str) // integer or floating point
{const char *nptr=str.c_str(),*last_character=nptr+str.size();char *endptr;
if(strtod(nptr,&endptr)){}return endptr==last_character;}

int Compare_Strings(const std::string &str1,const std::string &str2,bool case_sensitive)
{if(case_sensitive) return str1.compare(str2); // use the built-in comparison
else{std::string::const_iterator iter1,iter2;
    for(iter1=str1.begin(),iter2=str2.begin();iter1!=str1.end()&&iter2!=str2.end();iter1++,iter2++)
        if(::toupper(*iter1) != ::toupper(*iter2)) return ((::toupper(*iter1) < ::toupper(*iter2))?-1:1);
    std::string::size_type length1=str1.length(),length2=str2.length();
    return (length1==length2)?0:((length1<length2)?-1:1);}}

std::string toupper(const std::string& str)
{std::string str_copy=str;for(std::string::iterator iter=str_copy.begin();iter!=str_copy.end();iter++) *iter=::toupper(*iter);return str_copy;}

bool Parse_Integer_Range(const std::string& str,ARRAY<int>& integer_list)
{const char* c_str=str.c_str();
const char* endptr=c_str;
int start_val=(int)strtol(c_str,&const_cast<char*&>(endptr),10);
if(endptr==c_str) return false;
while(isspace(*endptr)) endptr++; // skip whitespace
if(*endptr=='\0'){integer_list.Append(start_val);return true;} // single value
else if(*endptr=='-'){
    const char* nextptr=endptr+1;
    int end_val=(int)strtol(nextptr,&const_cast<char*&>(endptr),10);
    if(endptr==nextptr) return false;
    while(isspace(*endptr))endptr++; // skip whitespace
    int step;
    if(*endptr=='\0')
        step=1;
    else if(*endptr==':'){
        nextptr=endptr+1;
        step=(int)strtol(nextptr,&const_cast<char*&>(endptr),10);
        if(endptr==nextptr) return false;
        while(isspace(*endptr))endptr++; // skip whitespace
        if(*endptr) return false;}
    else return false;
    for(int i=start_val;i<=end_val;i+=step) integer_list.Append(i);
    return true;}
else return false;}

// integer list format: [range],[range],[range],... where each [range] is either a single <number> or <number>-<number>, e.g. "1-3,4,7,10-11"
bool Parse_Integer_List(const std::string& str,ARRAY<int>& integer_list)
{integer_list.Remove_All();
std::string remaining_string=str;
while(!str.empty()){
    std::string token;
    std::string::size_type comma_pos=remaining_string.find(",");
    if(comma_pos!=std::string::npos){token=remaining_string.substr(0,comma_pos);remaining_string=remaining_string.substr(comma_pos+1);}
    else{token=remaining_string;remaining_string="";}
    if(!Parse_Integer_Range(token,integer_list))return false;}
return true;}

std::string Join(const std::string& separator,const ARRAY<std::string>& tokens)
{
    if(tokens.m==0) return "";
    std::string str=tokens(1);
    for(int i=2;i<=tokens.m;i++) str+=separator+tokens(i);
    return str;
}

void Split(const std::string& str,const std::string& separator,ARRAY<std::string>& tokens)
{
    std::string remaining_string=str;std::string::size_type separator_length=separator.size();
    while(!remaining_string.empty()){
        std::string token;
        std::string::size_type separator_position=remaining_string.find(separator);
        if(separator_position!=std::string::npos){
            tokens.Append(remaining_string.substr(0,separator_position));
            remaining_string=remaining_string.substr(separator_position+separator_length);}
        else{tokens.Append(remaining_string);remaining_string="";}}
}

bool Ends_With(const std::string& input,const std::string& test)
{
    int input_index=input.length(),test_index=test.length();
    for(;test_index>0 && input_index>0;) if(input[--input_index]!=test[--test_index]) return false;
    return test_index==0;
}

bool IEnds_With(const std::string& input,const std::string& test)
{
    int input_index=input.length(),test_index=test.length();
    for(;test_index>0 && input_index>0;) if(std::toupper(input[--input_index])!=std::toupper(test[--test_index])) return false;
    return test_index==0;
}

std::string string_sprintf(const char *format_string,...) // Assumes a max string length of 2048, since Windows doesn't support safety
{
    static char tmp[2048];
    va_list marker;
    va_start(marker,format_string);
    vsprintf(tmp,format_string,marker);
    va_end(marker);
    return tmp;
}

//#####################################################################
}
}
