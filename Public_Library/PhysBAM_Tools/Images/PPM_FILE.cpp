//#####################################################################
// Copyright 2002-2005, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Images/PPM_FILE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Read
//#####################################################################
template<class T> void PPM_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("PPM_FILE::Read not implemented");
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> void PPM_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("PPM_FILE::Read not implemented");
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<int d> void PPM_FILE<T>::
Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{  
    assert(image.domain.min_corner.x==1 && image.domain.min_corner.y==1);
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filename,true,false); // no compression
    *output<<"P6"<<std::endl;
    *output<<"# Generated using PPM_FILE::Write"<<std::endl;
    *output<<image.counts.x<<" "<<image.counts.y<<std::endl;
    *output<<255<<std::endl;
    for(int j=image.counts.y;j>=1;j--)for(int i=1;i<=image.counts.x;i++){VECTOR<unsigned char,d> pixel=IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,j));Write_Binary<T>(*output,pixel[1],pixel[2],pixel[3]);}
    delete output;
}
//#####################################################################
// Function Is_Supported
//#####################################################################
template<class T> bool PPM_FILE<T>::
Is_Supported()
{
    return true;
}
//#####################################################################
template class PPM_FILE<float>;
template void PPM_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,3> ,VECTOR<int,2> >&);
template void PPM_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,4> ,VECTOR<int,2> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PPM_FILE<double>;
template void PPM_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,3> ,VECTOR<int,2> >&);
template void PPM_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,4> ,VECTOR<int,2> >&);
#endif
#endif
