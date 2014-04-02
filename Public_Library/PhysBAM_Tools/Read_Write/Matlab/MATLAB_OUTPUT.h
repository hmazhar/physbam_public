//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATLAB_OUTPUT 
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __MATLAB_OUTPUT__
#define __MATLAB_OUTPUT__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <fstream>
namespace PhysBAM{
template<class TV> class GRID;
template<class T,int d> class VECTOR;
class MATLAB_OUTPUT
{
public: 
    int little_endian; // (1) little endian, (0) big endian;
    
    MATLAB_OUTPUT();
    ~MATLAB_OUTPUT();

    void Check_Endian();
    template<class T> void Convert_Bytes(T& data) const;
    template<class T> void Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,1> >& grid,const int stepnumber);
    template<class T> void Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,1> >& x,const int stepnumber);
    template<class T> void Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,2> >& grid,const int stepnumber);
    template<class T> void Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,2> >& x,const ARRAY<T,VECTOR<int,2> >& y,const int stepnumber);
    template<class T> void Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,3> >& grid,const int stepnumber);
    template<class T> void Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,3> >& x,const ARRAY<T,VECTOR<int,3> >& y,const ARRAY<T,VECTOR<int,3> >& z,const int stepnumber);
    template<class T> void Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,1> >& output,const int stepnumber);
    template<class T> void Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,2> >& output,const int stepnumber);
    template<class T> void Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,3> >& output,const int stepnumber);
    template<class T,int d> void Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,d> >& X,const int stepnumber);
};
}
#endif
#endif
