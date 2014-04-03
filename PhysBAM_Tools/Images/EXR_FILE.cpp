//#####################################################################
// Copyright 2004-2006, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXR_FILE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Images/EXR_FILE.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
using namespace PhysBAM;

#ifdef USE_OPENEXR

#include <ImathBox.h>
#include <ImfRgbaFile.h>
using namespace Imf;

//#####################################################################
// Function Write
//#####################################################################
template<class T> template<int d> void EXR_FILE<T>::
Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{
    assert(image.domain.min_corner.x==1&&image.domain.min_corner.y==1&&image.length==1);
    Rgba* pixels=new Rgba[image.counts.Product()];
    int t=0;
    if(d==3) for(int i=1;i<=image.counts.x;i++) for(int j=image.counts.y;j>=1;j--){pixels[t].r=(float)image(i,j)[1];pixels[t].g=(float)image(i,j)[2];pixels[t].b=(float)image(i,j)[3];pixels[t].a=(float)1;t++;}
    else if (d==4) for(int i=1;i<=image.counts.x;i++) for(int j=image.counts.y;j>=1;j--){pixels[t].r=(float)image(i,j)[1];pixels[t].g=(float)image(i,j)[2];pixels[t].b=(float)image(i,j)[3];pixels[t].a=(float)image(i,j)[4];t++;}
    try{RgbaOutputFile file(filename.c_str(),image.counts.x,image.counts.y,WRITE_RGBA);file.setFrameBuffer(pixels,image.counts.y,1);file.writePixels(image.counts.y);}
    catch(...){LOG::cerr<<"Cannot write exr image to "<<filename<<std::endl;PHYSBAM_FATAL_ERROR();}
    delete[] pixels;
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> void EXR_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image)
{
    try{
        RgbaInputFile file(filename.c_str());Imath::Box2i data_window=file.dataWindow();
        int width=data_window.max.x-data_window.min.x+1,height=data_window.max.y-data_window.min.y+1;
        Rgba* pixels=new Rgba[height*width];
        file.setFrameBuffer(pixels, 1, width);
        file.readPixels(data_window.min.y,data_window.max.y);
        image.Resize(1,width,1,height);
        int t=0;
        for(int j=image.counts.y;j>=1;j--) for(int i=1;i<=image.counts.x;i++){image(i,j)=VECTOR<T,3>((T)pixels[t].r,(T)pixels[t].g,(T)pixels[t].b);t++;}
        delete[] pixels;}
    catch(...){LOG::cerr<<"Cannot read exr image from "<<filename<<std::endl;PHYSBAM_FATAL_ERROR();}
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> void EXR_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("EXR_FILE: cannot read alpha");
}
//#####################################################################
// Function Is_Supported
//#####################################################################
template<class T> bool EXR_FILE<T>::
Is_Supported()
{
    return true;
}

#else
//#####################################################################
// Function Read
//#####################################################################
template<class T> void EXR_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("Not compiled with USE_OPENEXR.  Cannot read exr image.");
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> void EXR_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("Not compiled with USE_OPENEXR.  Cannot read exr image.");
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<int d> void EXR_FILE<T>::
Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("Not compiled with USE_OPENEXR.  Cannot write exr image.");
}
//#####################################################################
// Function Is_Supported
//#####################################################################
template<class T> bool EXR_FILE<T>::
Is_Supported()
{
    return false;
}
//#####################################################################
#endif
template class EXR_FILE<float>;
template void EXR_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,3> ,VECTOR<int,2> >&);
template void EXR_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,4> ,VECTOR<int,2> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EXR_FILE<double>;
template void EXR_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,3> ,VECTOR<int,2> >&);
template void EXR_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,4> ,VECTOR<int,2> >&);
#endif
#endif
