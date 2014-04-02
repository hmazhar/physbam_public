//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Read_Write/Matlab/MATLAB_OUTPUT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
MATLAB_OUTPUT::
MATLAB_OUTPUT()
{
    Check_Endian();
}
//#####################################################################
// Destructor
//#####################################################################
MATLAB_OUTPUT::
~MATLAB_OUTPUT()
{}
//#####################################################################
// Function Check_Endian
//#####################################################################
void MATLAB_OUTPUT::
Check_Endian()
{
    unsigned short one=0x0001;
    if(*((char*)&one)) little_endian=1; // little endian (e.g. PC), 01 comes before 00
    else little_endian=0;                    // big endian (e.g. SGI), 00 comes before 01
}
//#####################################################################
// Function Convert_Bytes
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Convert_Bytes(T& data) const
{
    union {T data;char raw[sizeof(T)];} input,output;
    input.data=data;
    for(int k=1;k<=(int)sizeof(T);k++) output.raw[k-1]=input.raw[sizeof(T)-k];
    data=output.data;
}
//#####################################################################
// Function Write_Header_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,1> >& grid,const int stepnumber)
{
    int m=grid.counts.x;
    ARRAY<T,VECTOR<int,1> > x(1,m);for(int i=1;i<=m;i++) x(i)=grid.Axis_X(i,1);
    Write_Header_File(file_name,x,stepnumber);
}
//#####################################################################
// Function Write_Header_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,1> >& x,const int stepnumber)
{
    int m=x.domain.max_corner.x,data_int;double data_double;
    std::ofstream Matlab_Output;
    Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    data_int=m;if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    for(int i=1;i<=m;i++){data_double=x(i);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Header_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,2> >& grid,const int stepnumber)
{
    int m=grid.Counts().x,n=grid.Counts().y;
    ARRAY<T,VECTOR<int,2> > x(1,m,1,n),y(1,m,1,n);for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){x(i,j)=grid.X(i,j).x;y(i,j)=grid.X(i,j).y;}
    Write_Header_File(file_name,x,y,stepnumber);
}
//#####################################################################
// Function Write_Header_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,2> >& x,const ARRAY<T,VECTOR<int,2> >& y,const int stepnumber)
{
    int m=x.domain.max_corner.x,n=x.domain.max_corner.y,data_int;double data_double;
    std::ofstream Matlab_Output;
    Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    data_int=m;if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    data_int=n;if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    int i;for(i=1;i<=m;i++) for(int j=1;j<=n;j++){data_double=x(i,j);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    for(i=1;i<=m;i++) for(int j=1;j<=n;j++){data_double=y(i,j);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Header_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,3> >& grid,const int stepnumber)
{
    int m=grid.Counts().x,n=grid.Counts().y,mn=grid.Counts().z;
    ARRAY<T,VECTOR<int,3> > x(1,m,1,n,1,mn),y(1,m,1,n,1,mn),z(1,m,1,n,1,mn);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){x(i,j,ij)=grid.X(i,j,ij).x;y(i,j,ij)=grid.X(i,j,ij).y;z(i,j,ij)=grid.X(i,j,ij).z;}
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){x(i,j,ij)=grid.X(i,j,ij).x;y(i,j,ij)=grid.X(i,j,ij).y;z(i,j,ij)=grid.X(i,j,ij).z;}
    Write_Header_File(file_name,x,y,z,stepnumber);
}
//#####################################################################
// Function Write_Header_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,3> >& x,const ARRAY<T,VECTOR<int,3> >& y,const ARRAY<T,VECTOR<int,3> >& z,const int stepnumber)
{
    int m=x.domain.max_corner.x,n=x.domain.max_corner.y,mn=x.domain.max_corner.z,data_int;double data_double;
    std::ofstream Matlab_Output;Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    data_int=m;if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    data_int=n;if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    data_int=mn;if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    int i;for(i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){
        data_double=x(i,j,ij);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    for(i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){
        data_double=y(i,j,ij);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    for(i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){
        data_double=z(i,j,ij);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,1> >& output,const int stepnumber)
{
    int m=output.domain.max_corner.x;double data_double;
    std::ofstream Matlab_Output;Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    for(int i=1;i<=m;i++){data_double=output(i);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,2> >& output,const int stepnumber)
{
    int m=output.domain.max_corner.x,n=output.domain.max_corner.y;double data_double;
    std::ofstream Matlab_Output;Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){data_double=output(i,j);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,3> >& output,const int stepnumber)
{
    int m=output.domain.max_corner.x,n=output.domain.max_corner.y,mn=output.domain.max_corner.z;double data_double;
    std::ofstream Matlab_Output;Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){
        data_double=output(i,j,ij);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T,int d> void MATLAB_OUTPUT::
Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,d> >& X,const int stepnumber)
{
    std::ofstream Matlab_Output;
    Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    int data_int=X.Size();if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    for(int a=1;a<=d;a++) for(int k=1;k<=X.Size();k++){
        double data_double=X(k)[a];if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template void MATLAB_OUTPUT::Convert_Bytes(T& data) const; \
    template void MATLAB_OUTPUT::Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,1> >& grid,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,1> >& x,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,2> >& grid,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,2> >& x,const ARRAY<T,VECTOR<int,2> >& y,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,3> >& grid,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,3> >& x,const ARRAY<T,VECTOR<int,3> >& y,const ARRAY<T,VECTOR<int,3> >& z,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,1> >& output,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,2> >& output,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,3> >& output,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,1> >& X,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,2> >& X,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,3> >& X,const int stepnumber);

INSTANTIATION_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double);
#endif
