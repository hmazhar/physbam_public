//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Read_Write/Matlab/GNUPLOT_OUTPUT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
GNUPLOT_OUTPUT::
GNUPLOT_OUTPUT()
{}
//#####################################################################
// Destructor
//#####################################################################
GNUPLOT_OUTPUT::
~GNUPLOT_OUTPUT()
{}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T> void GNUPLOT_OUTPUT::
Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,1> >& grid,const ARRAY<T,VECTOR<int,1> >& output,const int stepnumber)
{
    int m=output.domain.max_corner.x;
    std::ofstream Matlab_Output;Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str());
    for(int i=1;i<=m;i++){VECTOR<T,1> position=grid.X(i);Matlab_Output<<position.x<<"\t"<<output(i)<<std::endl;}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T> void GNUPLOT_OUTPUT::
Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,2> >& grid,const ARRAY<T,VECTOR<int,2> >& output,const int stepnumber)
{
    int m=output.domain.max_corner.x,n=output.domain.max_corner.y;
    std::ofstream Matlab_Output;Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str());
    for(int i=1;i<=m;i++){for(int j=1;j<=n;j++){
            VECTOR<T,2> position=grid.X(i,j);Matlab_Output<<position.x<<"\t"<<position.y<<"\t"<<output(i,j)<<std::endl;}
        Matlab_Output<<std::endl;}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T> void GNUPLOT_OUTPUT::
Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,3> >& grid,const ARRAY<T,VECTOR<int,3> >& output,const int stepnumber)
{
    int m=output.domain.max_corner.x,n=output.domain.max_corner.y,mn=output.domain.max_corner.z;
    std::ofstream Matlab_Output;Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str());
    for(int i=1;i<=m;i++){for(int j=1;j<=n;j++){for(int ij=1;ij<=mn;ij++){
            VECTOR<T,3> position=grid.X(i,j,ij);Matlab_Output<<position.x<<"\t"<<position.y<<"\t"<<position.z<<"\t"<<output(i,j,ij)<<std::endl;}
            Matlab_Output<<std::endl;}
        Matlab_Output<<std::endl;}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
#ifdef COMPILE_WITH_BINTREE_SUPPORT
template<class T> void GNUPLOT_OUTPUT::
Write_Output_File(const std::string& file_name,const BINTREE_GRID<T>& grid,const ARRAY<T>& output,const int stepnumber)
{
    ARRAY<PAIR<T,int> > list_to_sort;
    for(DYADIC_GRID_ITERATOR_CELL<BINTREE_GRID<T> > iterator(grid,0);iterator.Valid();iterator.Next()){
        list_to_sort.Append(PAIR<T,int>(iterator.Location().x,iterator.Cell_Index()));}
    Sort(list_to_sort);

    std::ofstream Matlab_Output;Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str());
    for(int i=1;i<=list_to_sort.Size();++i) Matlab_Output<<list_to_sort(i).x<<"\t"<<output(list_to_sort(i).y)<<std::endl;
    Matlab_Output.close();
}
#endif
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T,int d> void GNUPLOT_OUTPUT::
Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,d> >& X,const int stepnumber)
{
    std::ofstream Matlab_Output;
    Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str());
    for(int i=1;i<=X.Size();++i){for(int a=1;a<=d;++a){
        Matlab_Output<<X(i)[a];if(a!=d) Matlab_Output<<"\t";}
        Matlab_Output<<std::endl;}
    Matlab_Output.close();
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template void GNUPLOT_OUTPUT::Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,1> >& grid,const ARRAY<T,VECTOR<int,1> >& output,const int stepnumber); \
    template void GNUPLOT_OUTPUT::Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,2> >& grid,const ARRAY<T,VECTOR<int,2> >& output,const int stepnumber); \
    template void GNUPLOT_OUTPUT::Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,3> >& grid,const ARRAY<T,VECTOR<int,3> >& output,const int stepnumber); \
    template void GNUPLOT_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,1> >& X,const int stepnumber); \
    template void GNUPLOT_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,2> >& X,const int stepnumber); \
    template void GNUPLOT_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,3> >& X,const int stepnumber);

INSTANTIATION_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double);
#endif

#ifdef COMPILE_WITH_BINTREE_SUPPORT
template void GNUPLOT_OUTPUT::Write_Output_File(const std::string&,const BINTREE_GRID<float>&,const ARRAY<float>&,const int);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void GNUPLOT_OUTPUT::Write_Output_File(const std::string&,const BINTREE_GRID<double>&,const ARRAY<double>&,const int);
#endif
#endif
