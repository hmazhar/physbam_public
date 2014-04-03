//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMAGE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Images/BMP_FILE.h>
#include <PhysBAM_Tools/Images/EXR_FILE.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Images/JPG_FILE.h>
#include <PhysBAM_Tools/Images/PNG_FILE.h>
#include <PhysBAM_Tools/Images/PPM_FILE.h>
#include <PhysBAM_Tools/Images/RGB_FILE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Function Read
//#####################################################################
template<class T> template<int d> void IMAGE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{
    std::string extension=FILE_UTILITIES::Get_File_Extension(filename);
    if(extension=="bmp") BMP_FILE<T>::Read(filename,image);
    else if(extension=="jpg") JPG_FILE<T>::Read(filename,image);
    else if(extension=="ppm") PPM_FILE<T>::Read(filename,image);
    else if(extension=="rgb") RGB_FILE<T>::Read(filename,image);
    else if(extension=="png") PNG_FILE<T>::Read(filename,image);
    else if(extension=="exr") EXR_FILE<T>::Read(filename,image);
    else if(extension=="pbi") FILE_UTILITIES::Read_From_File<float>(filename,image);
    else PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unknown image file extension  from filename '%s' extension '%s'",filename.c_str(),extension.c_str()));
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<int d> void IMAGE<T>::
Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image,const T gamma,const T dither_amplitude)
{
    RANDOM_NUMBERS<T> random;random.Set_Seed(324032); // want uniform seed so noise perturbation pattern is temporally coherent
    ARRAY<VECTOR<T,d> ,VECTOR<int,2> > *corrected_image=0; // may need to shift image to be (1,1) based, or to gamma correct it
    if(gamma!=1 || image.domain.min_corner.x!=1 || image.domain.min_corner.y!=1 || dither_amplitude>0){
        corrected_image=new ARRAY<VECTOR<T,d> ,VECTOR<int,2> >(1,image.counts.x,1,image.counts.y,false);
        ARRAY<VECTOR<T,d> ,VECTOR<int,2> >::Shifted_Get(*corrected_image,image,VECTOR<int,2>(image.domain.min_corner.x-1,image.domain.min_corner.y-1));
        T one_over_gamma=1/gamma;
        for(int t=1;t<=corrected_image->array.Size();t++){
            VECTOR<T,d> color=corrected_image->array(t);
            for(int channel=1;channel<=d;channel++) corrected_image->array(t)[channel]=pow(color[channel],one_over_gamma);
            if(dither_amplitude){
                VECTOR<T,d> pixel_values((T)256*corrected_image->array(t));
                VECTOR<int,d> floored_values;
                for(int channel=1;channel<=d;channel++) floored_values[channel]=(int)pixel_values[channel];
                VECTOR<T,d> random_stuff=random.Get_Uniform_Vector(VECTOR<T,d>(),VECTOR<T,d>::All_Ones_Vector());
                VECTOR<T,d> normalized_values=pixel_values-VECTOR<T,d>(floored_values);
                for(int k=1;k<=d;k++)
                    if(random_stuff(k)>normalized_values(k)) corrected_image->array(t)[k]=(floored_values[k]+(T).5001)/(T)256; // use normal quantized floor
                    else corrected_image->array(t)[k]=(floored_values[k]+(T)1.5001)/(T)256;}}} // jump to next value
    const ARRAY<VECTOR<T,d> ,VECTOR<int,2> > &image_to_write=corrected_image?*corrected_image:image;
    std::string extension=FILE_UTILITIES::Get_File_Extension(filename);
    if(extension=="bmp") BMP_FILE<T>::Write(filename,image_to_write);
    else if(extension=="jpg") JPG_FILE<T>::Write(filename,image_to_write);
    else if(extension=="ppm") PPM_FILE<T>::Write(filename,image_to_write);
    else if(extension=="png") PNG_FILE<T>::Write(filename,image_to_write);
    else if(extension=="rgb") RGB_FILE<T>::Write(filename,image_to_write);
    else if(extension=="exr") EXR_FILE<T>::Write(filename,image_to_write);
    else if(extension=="pbi") FILE_UTILITIES::Write_To_File<float>(filename,image_to_write);
    else PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unknown image file extension from filename '%s' extension '%s'",filename.c_str(),extension.c_str()));
    delete corrected_image;
}
//#####################################################################
// Function Is_Supported
//#####################################################################
template<class T> bool IMAGE<T>::
Is_Supported(const std::string& filename)
{
    std::string extension=FILE_UTILITIES::Get_File_Extension(filename);
    if(extension=="bmp") return BMP_FILE<T>::Is_Supported();
    else if(extension=="jpg") return JPG_FILE<T>::Is_Supported();
    else if(extension=="ppm") return PPM_FILE<T>::Is_Supported();
    else if(extension=="png") return PNG_FILE<T>::Is_Supported();
    else if(extension=="rgb") return RGB_FILE<T>::Is_Supported();
    else if(extension=="exr") return EXR_FILE<T>::Is_Supported();
    else if(extension=="pbi") return true;
    else return false;
}
//#####################################################################
template class IMAGE<float>;
template void IMAGE<float>::Read(const std::string&,ARRAY<VECTOR<float,3> ,VECTOR<int,2> >&);
template void IMAGE<float>::Read(const std::string&,ARRAY<VECTOR<float,4> ,VECTOR<int,2> >&);
template void IMAGE<float>::Write(const std::string&,const ARRAY<VECTOR<float,3> ,VECTOR<int,2> >&,const float,const float);
template void IMAGE<float>::Write(const std::string&,const ARRAY<VECTOR<float,4> ,VECTOR<int,2> >&,const float,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class IMAGE<double>;
template void IMAGE<double>::Read(const std::string&,ARRAY<VECTOR<double,3> ,VECTOR<int,2> >&);
template void IMAGE<double>::Read(const std::string&,ARRAY<VECTOR<double,4> ,VECTOR<int,2> >&);
template void IMAGE<double>::Write(const std::string&,const ARRAY<VECTOR<double,3> ,VECTOR<int,2> >&,const double,const double);
template void IMAGE<double>::Write(const std::string&,const ARRAY<VECTOR<double,4> ,VECTOR<int,2> >&,const double,const double);
#endif
#endif
