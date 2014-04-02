//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_BMP_HEADER
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_BMP_HEADER__
#define __READ_WRITE_BMP_HEADER__

#include <PhysBAM_Tools/Images/BMP_HEADER.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW>
class Read_Write<BMP_HEADER,RW>
{
public:
    static void Read(std::istream& input,BMP_HEADER& object)
    {Read_Binary<RW>(input,object.file_type[0],object.file_type[1],object.file_size,object.reserved1,object.reserved2,object.offset,object.info_header_size);
    Read_Binary<RW>(input,object.w,object.h,object.number_of_bitplanes,object.bits_per_pixel,object.type_of_compression);
    Read_Binary<RW>(input,object.bitmap_size,object.x_pixels_per_meter,object.y_pixels_per_meter,object.number_of_colors,object.number_of_important_colors);

    //check validity
    if(object.file_type[0]!='B' || object.file_type[1]!='M') PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Illegal file type: %c%c",object.file_type[0],object.file_type[1]));
    if(object.info_header_size != 40) LOG::cerr<<"Warning: weird info_header_size: "<<object.info_header_size<<" (expected 40)"<<std::endl;
    if(object.number_of_bitplanes != 1) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Illegal number of bitplanes: %d (must be 1)",object.number_of_bitplanes));
    if(object.bits_per_pixel != 24) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Number of bits per pixel: %d (the only supported number is 24)",object.bits_per_pixel));
    if(object.type_of_compression != 0) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Type of compression %d (the only supported type is 0)",object.type_of_compression));}

    static void Write(std::ostream& output,const BMP_HEADER& object)
    {Write_Binary<RW>(output,object.file_type[0],object.file_type[1],object.file_size,object.reserved1,object.reserved2,object.offset,object.info_header_size);
    Write_Binary<RW>(output,object.w,object.h,object.number_of_bitplanes,object.bits_per_pixel,object.type_of_compression);
    Write_Binary<RW>(output,object.bitmap_size,object.x_pixels_per_meter,object.y_pixels_per_meter,object.number_of_colors,object.number_of_important_colors);}
};
template<class T>
inline std::ostream& operator<<(std::ostream& output,const BMP_HEADER& header)
{output<<"BMP_HEADER:"<<std::endl
    <<"file_type: "<<header.file_type[0]<<" "<<header.file_type[1]<<std::endl
    <<"file_size  "<<int(header.file_size)<<std::endl
    <<"reserved1  "<<int(header.reserved1)<<std::endl
    <<"reserved2  "<<int(header.reserved2)<<std::endl
    <<"offset     "<<int(header.offset)<<std::endl
    <<"info_header_size    "<<header.info_header_size<<std::endl
    <<"w "<<header.w<<", h "<<header.h<<std::endl
    <<"number_of_bitplanes "<<int(header.number_of_bitplanes)<<std::endl
    <<"bits_per_pixel      "<<int(header.bits_per_pixel)<<std::endl
    <<"type_of_compression "<<header.type_of_compression<<std::endl
    <<"bitmap_size         "<<header.bitmap_size<<std::endl
    <<"x_pixels_per_meter  "<<header.x_pixels_per_meter<<std::endl
    <<"y_pixels_per_meter  "<<header.y_pixels_per_meter<<std::endl
    <<"number_of_colors    "<<header.number_of_colors<<std::endl
    <<"number_of_important_colors "<<header.number_of_important_colors<<std::endl;
return output;}
}
#endif
#endif
