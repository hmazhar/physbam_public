//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RGB_HEADER
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RGB_HEADER__
#define __READ_WRITE_RGB_HEADER__

#include <PhysBAM_Tools/Images/RGB_HEADER.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Utilities/EXCEPTIONS.h>
namespace PhysBAM{

template<class RW>
class Read_Write<RGB_HEADER,RW>
{
public:
    static void Read(std::istream& input,RGB_HEADER& object)
    {Read_Binary<RW>(input,object.magic_number,object.compression,object.bytes_per_channel,object.dimensions,object.width,object.height);
    Read_Binary<RW>(input,object.channels,object.pixel_minimum_value,object.pixel_maximum_value,object.dummy);
    Read_Binary_Array<RW>(input,object.name,80);Read_Binary<RW>(input,object.colormap);Read_Binary_Array<RW>(input,object.dummy2,404);
    object.Swap_Endian();
    
    //check validity
    if(object.magic_number!=474) throw READ_ERROR(STRING_UTILITIES::string_sprintf("illegal file type %d, should be 474",object.magic_number));
    if(object.channels != 3) throw READ_ERROR(STRING_UTILITIES::string_sprintf("illegal number of channels: %d (must be 3)",object.channels));
    if(object.bytes_per_channel != 1) throw READ_ERROR(STRING_UTILITIES::string_sprintf("only 1 bytes per channel supported (you gave %d)",object.bytes_per_channel));
    if(object.dimensions != 2 && object.dimensions != 3) throw READ_ERROR(STRING_UTILITIES::string_sprintf("dimensions field must be 2 or 3 (you gave %d)",object.dimensions));}

    static void Write(std::ostream& output,const RGB_HEADER& object)
    {RGB_HEADER swapped=object;swapped.Swap_Endian();
    Write_Binary<RW>(output,swapped.magic_number,swapped.compression,swapped.bytes_per_channel,swapped.dimensions,swapped.width,swapped.height,swapped.channels);
    Write_Binary<RW>(output,swapped.pixel_minimum_value,swapped.pixel_maximum_value,swapped.dummy);
    Write_Binary_Array<RW>(output,object.name,80);Write_Binary<RW>(output,object.colormap);Write_Binary_Array<RW>(output,object.dummy2,404);}
};
}
#endif
#endif
