//#####################################################################
// Copyright 2003s-2005, Geoffrey Irving, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class JPG_FILE  
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Images/JPG_FILE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
using namespace PhysBAM;

#ifdef USE_LIBJPEG

#include <cstdio>
extern "C"{
#include <jpeglib.h>
}
//#####################################################################
// Function Read_Error
//#####################################################################
static void Read_Error(j_common_ptr cinfo)
{
    throw READ_ERROR("JPG_FILE:: Can't read image");
}
template<class T> void JPG_FILE<T>::
//#####################################################################
// Function Read
//#####################################################################
Read(const std::string& filename,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& image)
{
    throw READ_ERROR("JPG_FILE:: Cannot read alpha");
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> void JPG_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image)
{
    struct jpeg_decompress_struct cinfo;FILE * infile;int row_stride;struct jpeg_error_mgr error_manager;
    if(!(infile=fopen(filename.c_str(),"rb"))) throw READ_ERROR(STRING_UTILITIES::string_sprintf("JPG_FILE::Read: Can't open %s",filename.c_str()));
    cinfo.err=jpeg_std_error(&error_manager);error_manager.error_exit=Read_Error;
    jpeg_create_decompress(&cinfo);jpeg_stdio_src(&cinfo,infile);jpeg_read_header(&cinfo,TRUE);jpeg_start_decompress(&cinfo);

    row_stride=cinfo.output_width*cinfo.output_components;
    JSAMPLE* row=new unsigned char[row_stride];JSAMPROW row_pointer[]={row};
    LOG::cout<<"reading "<<filename<<": "<<row_stride/3<<" x "<<cinfo.output_height<<std::endl;

    image.Resize(1,cinfo.output_width,1,cinfo.output_height,false,false);
    while(cinfo.output_scanline<cinfo.output_height){
        jpeg_read_scanlines(&cinfo,row_pointer,1);int index=0;
        for(int i=1;i<=image.counts.x;i++){
            unsigned char r=row[index++],g=row[index++],b=row[index++];
            image(i,image.counts.y-cinfo.output_scanline+1)=IMAGE<T>::Byte_Color_To_Scalar_Color(VECTOR<unsigned char,3>(r,g,b));}}
    jpeg_finish_decompress(&cinfo);jpeg_destroy_decompress(&cinfo);delete[] row;

    fclose(infile);
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<int d> void JPG_FILE<T>::
Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{  
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE* outfile; // target file

    cinfo.err=jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    if(!(outfile = fopen(filename.c_str(), "wb")))
        PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("JPG_FILE::Write: Can't open %s",filename.c_str()));
    jpeg_stdio_dest(&cinfo,outfile);
    cinfo.image_width=image.counts.x;cinfo.image_height=image.counts.y;cinfo.input_components=3;
    cinfo.in_color_space=JCS_RGB; // colorspace of input image
    jpeg_set_defaults(&cinfo);jpeg_set_quality(&cinfo,95,TRUE); // limit to baseline-JPEG values
    jpeg_start_compress(&cinfo,TRUE);

    int row_stride=cinfo.image_width*3; // JSAMPLEs per row in image_buffer
    JSAMPLE* row=new unsigned char[row_stride];JSAMPROW row_pointer[]={row};
    while(cinfo.next_scanline < cinfo.image_height){
        int index=0;for(int i=1;i<=image.counts.x;i++){VECTOR<unsigned char,d> pixel=IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,image.counts.y-cinfo.next_scanline));row[index++]=pixel[1];row[index++]=pixel[2];row[index++]=pixel[3];} // copy row
        jpeg_write_scanlines(&cinfo,row_pointer,1);}
    delete[] row;
    jpeg_finish_compress(&cinfo);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
}
//#####################################################################
// Function Is_Supported
//#####################################################################
template<class T> bool JPG_FILE<T>::
Is_Supported()
{
    return true;
}

#else

//#####################################################################
// Function Read
//#####################################################################
template<class T> void JPG_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("Not compiled with USE_LIBJPEG.  Cannot read jpeg image.");
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> void JPG_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("Not compiled with USE_LIBJPEG.  Cannot read jpeg image.");
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<int d> void JPG_FILE<T>::
Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{  
    PHYSBAM_FATAL_ERROR("Not compiled with USE_LIBJPEG.  Cannot write jpeg image.");
}
//#####################################################################
// Function Is_Supported
//#####################################################################
template<class T> bool JPG_FILE<T>::
Is_Supported()
{
    return false;
}
//#####################################################################
#endif

template class JPG_FILE<float>;
template void JPG_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,3> ,VECTOR<int,2> >&);
template void JPG_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,4> ,VECTOR<int,2> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class JPG_FILE<double>;
template void JPG_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,3> ,VECTOR<int,2> >&);
template void JPG_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,4> ,VECTOR<int,2> >&);
#endif
#endif
