//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace FILE_UTILITIES
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/ZIP.h>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <fstream>
#ifdef WIN32
#include <windows.h>
#endif
#if defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__)
#include <sys/stat.h>
#endif
namespace PhysBAM{
namespace FILE_UTILITIES{

static const char *file_extensions[]={"rgd","rgd2d","tri","phi","phi2d","oct","ply","ply2d","rgb","tri2d","curve","curve2d","tet","hex","box","phoneme",0};

//###################################################################
// Win32 Specific Function Definitions
//###################################################################
#if defined(WIN32)

bool Directory_Exists(const std::string& dirname)
{DWORD attr=GetFileAttributes(dirname.c_str());return((attr!=-1)&&(attr&FILE_ATTRIBUTE_DIRECTORY));}

bool Create_Directory(const std::string& dirname)
{if(!Directory_Exists(dirname)){
    LOG::cout<<"Creating directory using CreateDirectory...";CreateDirectory(dirname.c_str(),0);
    if(!Directory_Exists(dirname)){LOG::cerr<<"Failed!"<<std::endl;throw FILESYSTEM_ERROR("Create_Directory failed");return false;}
    LOG::cout<<"Successful!"<<std::endl;}
return true;}

std::string Real_Path(const std::string& path)//TODO: Implement this for windows
{PHYSBAM_NOT_IMPLEMENTED();}

int Compare_File_Times_Ignoring_Compression_Suffix(const std::string& filename1,const std::string& filename2)
{HANDLE handle1=CreateFile(filename1.c_str(),0,0,0,OPEN_EXISTING,0,0);
if(handle1==INVALID_HANDLE_VALUE){LOG::cerr<<"Compare_File_Times: can't CreateFile "<<filename1<<std::endl;return 0;}
HANDLE handle2=CreateFile(filename2.c_str(),0,0,0,OPEN_EXISTING,0,0);
if(handle2==INVALID_HANDLE_VALUE){LOG::cerr<<"Compare_File_Times: can't CreateFile "<<filename2<<std::endl;return 0;}
FILETIME time1,time2;
if(!GetFileTime(handle1,0,&time1,0)||!GetFileTime(handle2,0,&time2,0)){LOG::cerr<<"Compare_File_Times: error with GetFileTime"<<std::endl;return 0;}
CloseHandle(handle1);CloseHandle(handle2);return CompareFileTime(&time1,&time2);}

FILE* Temporary_File()
{return fopen(_tempnam(getenv("TEMP"),"pb"),"w");}

std::string Get_Working_Directory()
{return ".";}

//###################################################################
// Linux Specific Function Definitions
//###################################################################
#elif defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__)

bool Directory_Exists(const std::string& dirname)
{return std::ifstream(dirname.c_str())!=0;}

bool Create_Directory(const std::string& dirname)
{if(!Directory_Exists(dirname)){
    std::string command="mkdir -p "+dirname;
    LOG::cout<<"Creating directory using system(\""<<command<<"\")...";
    if(system(command.c_str()) || !Directory_Exists(dirname)){LOG::cerr<<"Failed!"<<std::endl;throw FILESYSTEM_ERROR("Create_Directory failed");return false;}
    LOG::cout<<"Successful!"<<std::endl;}
return true;}

std::string Real_Path(const std::string& path)
{char real_path[PATH_MAX];PHYSBAM_ASSERT(realpath(path.c_str(),real_path));return real_path;}

int Compare_File_Times_Ignoring_Compression_Suffix(const std::string& filename1,const std::string& filename2)
{struct stat stat1,stat2;
if(stat(filename1.c_str(),&stat1)!=0){LOG::cerr<<"Compare_File_Times: can't stat "<<filename1<<std::endl;}
if(stat(filename2.c_str(),&stat2)!=0){LOG::cerr<<"Compare_File_Times: can't stat "<<filename2<<std::endl;}
if(stat1.st_mtime<stat2.st_mtime)return -1;
else if(stat1.st_mtime>stat2.st_mtime)return 1;
else return 0;}

FILE* Temporary_File()
{return tmpfile();}

std::string Get_Working_Directory()
{ARRAY<char> buffer(128,false);
for(;;){
    if(getcwd(buffer.Get_Array_Pointer(),buffer.m-1))
        return std::string(buffer.Get_Array_Pointer());
    else if(errno==ERANGE){
        if(buffer.m>=4096) PHYSBAM_FATAL_ERROR("refusing to allocate more than 4k to return working directory");
        buffer.Resize(2*buffer.m,false,false);}
    else PHYSBAM_FATAL_ERROR();}}

//###################################################################
// Default (Unimplemented) Function Definitions
//###################################################################
#else

bool Directory_Exists(const std::string& dirname)
{return true;}  // always return true on unsupported platforms

std::string Real_Path(const std::string& path)
{PHYSBAM_NOT_IMPLEMENTED();}

int Compare_File_Times_Ignoring_Compression_Suffix(const std::string& filename1,const std::string& filename2)
{PHYSBAM_NOT_IMPLEMENTED();}

FILE* Temporary_File()
{PHYSBAM_NOT_IMPLEMENTED();}

std::string Get_Working_Directory()
{return ".";}

#endif
//###################################################################
// File open with compression capability
//###################################################################
std::istream* Safe_Open_Input(const std::string& filename,bool binary)
{
    bool compressed=File_Is_Compressed(filename);
    std::ios_base::openmode flags=std::ios::in;if(binary) flags|=std::ios::binary;
    std::string filename_compressed=compressed?filename:filename+".gz";
#ifndef COMPILE_WITHOUT_ZLIB_SUPPORT
    if(File_Exists(filename_compressed))
        return Gzip_In(filename_compressed,flags);
#else
    if(File_Exists(filename_compressed))
        PHYSBAM_FATAL_ERROR("Must compile with zlib to read compressed files");
#endif
    if(!compressed){
        std::istream* input=new std::ifstream(filename.c_str(),flags);if(*input) return input;
        filename_compressed=filename+"(.gz)";}
    LOG::cerr<<"Can't open "<<filename_compressed<<" for read "<<(binary?"(binary)":"")<<std::endl;
    throw FILESYSTEM_ERROR("Safe_Open_Input failed");
}

std::ostream* Safe_Open_Output(const std::string& filename,bool binary,bool write_compressed_if_possible)
{
    bool compressed=File_Is_Compressed(filename);
    std::ios_base::openmode flags=std::ios::out;if(binary) flags|=std::ios::binary; 
    if(!write_compressed_if_possible && File_Exists_Ignoring_Compression_Suffix(filename+".gz")){
        LOG::cerr<<"Refusing to write "<<filename<<" uncompressed when compressed version already exists\n";
        throw FILESYSTEM_ERROR("Safe_Open_Output failed (compressed already exists)");}
    if(compressed && !binary){LOG::cerr<<"Refusing to open compressed file "<<filename<<"in text mode\n";
        throw FILESYSTEM_ERROR("Safe_Open_Output failed (tried to make text compressed)");}
    std::string actual_filename;
#ifndef COMPILE_WITHOUT_ZLIB_SUPPORT
    if(binary && write_compressed_if_possible){actual_filename=compressed?filename:filename+".gz";
        if(File_Writable(actual_filename))
            return Gzip_Out(actual_filename,flags);
    }
    else{actual_filename=compressed?Strip_Compression_Suffix(filename):filename;
        std::ostream* output=new std::ofstream(actual_filename.c_str(),flags);if(*output) return output;}
#else
    actual_filename=compressed?Strip_Compression_Suffix(filename):filename;
    std::ostream* output=new std::ofstream(actual_filename.c_str(),flags);if(*output) return output;
#endif
    LOG::cerr<<"Can't open "<<actual_filename<<" for write "<<(binary?"(binary)":"")<<std::endl;
    throw FILESYSTEM_ERROR("Safe_Open_Input failed");
}
//###################################################################
// Function Compare_File_Times
//###################################################################
int Compare_File_Times(const std::string& filename1,const std::string& filename2)
{return Compare_File_Times_Ignoring_Compression_Suffix(Real_File(filename1),Real_File(filename2));}
//###################################################################
// Function File_Exists_Ignoring_Compression_Suffix
//###################################################################
bool File_Exists_Ignoring_Compression_Suffix(const std::string& filename)
{return std::ifstream(filename.c_str())!=0;}
//###################################################################
// Function File_Writable_Ignoring_Compression_Suffix
//###################################################################
bool File_Writable_Ignoring_Compression_Suffix(const std::string& filename)
{return std::ofstream(filename.c_str(),std::ios::out)!=0;} // TODO: make this not create the file
//###################################################################
// Function Directory_Writable
//###################################################################
bool Directory_Writable(const std::string& dirname) // TODO: make this nicer
{static const char* dummy_filename="_PHYSBAM_FILE_UTILITIES_DUMMY_";
std::string filename=dirname+"/"+dummy_filename;bool success=(std::ofstream(filename.c_str())!=0);remove(filename.c_str());return success;}
//###################################################################
// Function Find_First_Nonexistent_File_In_Sequence
//###################################################################
std::string Find_First_Nonexistent_File_In_Sequence(std::string filename_pattern,const int id_start,int* id_result)
{int id=id_start;while(File_Exists(STRING_UTILITIES::string_sprintf(filename_pattern.c_str(),id))) id++;
if(id_result) *id_result=id;return STRING_UTILITIES::string_sprintf(filename_pattern.c_str(),id);}
//###################################################################
// Function Find_First_Nonexistent_Directory_In_Sequence
//###################################################################
std::string Find_First_Nonexistent_Directory_In_Sequence(std::string directory_pattern,const int id_start,int* id_final)
{int id=id_start;while(Directory_Exists(STRING_UTILITIES::string_sprintf(directory_pattern.c_str(),id))) id++;
if(id_final) *id_final=id;return STRING_UTILITIES::string_sprintf(directory_pattern.c_str(),id);}
//###################################################################
// Function Make_First_Nonexistent_Directory_In_Sequence
//###################################################################
std::string Make_First_Nonexistent_Directory_In_Sequence(std::string directory_pattern,const int id_start,int* id_final)
{std::string output_directory=Find_First_Nonexistent_Directory_In_Sequence(directory_pattern,id_start,id_final);
Create_Directory(output_directory);return output_directory;}
//#####################################################################
std::string Get_File_Extension_Ignoring_Compression_Suffix(const std::string &filename)
{std::string::size_type lastperiod=filename.rfind('.'),lastslash=filename.rfind('/');
if(lastperiod!=std::string::npos&&(lastslash==std::string::npos||lastperiod>lastslash))return filename.substr(lastperiod+1);else return "";}

std::string Get_Basename_Ignoring_Compression_Suffix(const std::string& filename)
{std::string::size_type lastperiod=filename.rfind(".");std::string::size_type lastslash=filename.rfind("/");
if(lastperiod!=std::string::npos&&(lastslash==std::string::npos||lastperiod>lastslash))return filename.substr(0,lastperiod);else return filename;}

std::string Get_Short_Name_Ignoring_Compression_Suffix(const std::string& filename)
{size_t lastslash=filename.rfind('/');if(lastslash==std::string::npos)return filename;else return filename.substr(lastslash+1);}

bool File_Extension_Matches_Ignoring_Compression_Suffix(const std::string& filename,const std::string& ext,const bool case_sensitive)
{return !STRING_UTILITIES::Compare_Strings(Get_File_Extension_Ignoring_Compression_Suffix(filename),ext,case_sensitive);}

bool File_Is_Compressed(const std::string& filename)
{return File_Extension_Matches_Ignoring_Compression_Suffix(filename,"gz");}

bool File_Exists(const std::string& filename)
{return File_Exists_Ignoring_Compression_Suffix(filename) || (!File_Is_Compressed(filename) && File_Exists_Ignoring_Compression_Suffix(filename+".gz"));}

bool File_Writable(const std::string& filename)
{return File_Writable_Ignoring_Compression_Suffix(filename) || (!File_Is_Compressed(filename) && File_Writable_Ignoring_Compression_Suffix(filename+".gz"));}

std::string Strip_Compression_Suffix(const std::string& filename)
{if(File_Is_Compressed(filename)) return Get_Basename_Ignoring_Compression_Suffix(filename); else return filename;}

std::string Real_File(const std::string& filename)
{if(File_Exists_Ignoring_Compression_Suffix(filename))return filename;
else if(!File_Is_Compressed(filename) && File_Exists_Ignoring_Compression_Suffix(filename+".gz"))return filename+".gz";
else return "";}

bool File_Extension_Matches(const std::string& filename,const std::string& ext,const bool case_sensitive)
{return File_Extension_Matches_Ignoring_Compression_Suffix(Strip_Compression_Suffix(filename),ext,case_sensitive);}

FILE_TYPE Get_File_Type_Ignoring_Compression_Suffix(const std::string& filename)
{std::string ext=Get_File_Extension_Ignoring_Compression_Suffix(filename);
for(int i=0;file_extensions[i];i++)if(ext==file_extensions[i])return (FILE_TYPE)i;
return UNKNOWN_FILE;}

FILE_TYPE Get_File_Type(const std::string& filename)
{return Get_File_Type_Ignoring_Compression_Suffix(Strip_Compression_Suffix(filename));}

bool File_Type_Matches_Ignoring_Compression_Suffix(const std::string& filename,FILE_TYPE type)
{return File_Extension_Matches_Ignoring_Compression_Suffix(filename,file_extensions[(int)type]);}

bool File_Type_Matches(const std::string& filename,FILE_TYPE type)
{return File_Type_Matches_Ignoring_Compression_Suffix(Strip_Compression_Suffix(filename),type);}

std::string Get_File_Extension(const std::string &filename)
{return Get_File_Extension_Ignoring_Compression_Suffix(Strip_Compression_Suffix(filename));}

std::string Get_Basename(const std::string& filename)
{return Get_Basename_Ignoring_Compression_Suffix(Strip_Compression_Suffix(filename));}

std::string Get_Base_Directory_Name(const std::string& path)
{size_t lastslash=path.rfind('/');if(lastslash==std::string::npos)return ".";else return path.substr(0,lastslash);}

std::string Get_Short_Name(const std::string& filename)
{return Get_Short_Name_Ignoring_Compression_Suffix(Strip_Compression_Suffix(filename));}

bool Is_Animated(const std::string &filename)
{return filename.find("%d") != std::string::npos;}

std::string Get_Frame_Filename(const std::string &filename,int frame)
{return Is_Animated(filename)?STRING_UTILITIES::string_sprintf(filename.c_str(),frame):filename;}

bool Frame_File_Exists(const std::string &filename,int frame)
{return File_Exists(Get_Frame_Filename(filename,frame));}

void Ignore(std::istream& input,char c)
{
    while(isspace(input.peek())) input.get();
    if(input.peek()==c) input.get();
}

}
}
#endif
