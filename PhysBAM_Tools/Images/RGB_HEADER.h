//#####################################################################
// Copyright 2002-2005, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Struct RGB_HEADER
//#####################################################################
// A group of fields with the layout of a Windows rgb file header. Used by the class RGB_FILE.
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __RGB_HEADER__
#define __RGB_HEADER__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <cstring>
namespace PhysBAM{

struct RGB_HEADER
{
    unsigned short magic_number; // must be 474
    unsigned char compression; // 0 or 1
    unsigned char bytes_per_channel; // 1
    unsigned short dimensions; // 2 for a single 2D image
    unsigned short width,height;
    unsigned short channels; // 3 is RGB, 4 is RGBA
    unsigned int pixel_minimum_value; // typically 0
    unsigned int pixel_maximum_value; // typically 255
    unsigned int dummy; // ignored
    char name[80];
    unsigned int colormap;
    char dummy2[404];
    
    void Initialize(const int width_input,const int height_input)
    {magic_number=474;compression=0;bytes_per_channel=1;dimensions=2;width=width_input;height=height_input;channels=3;pixel_maximum_value=0;pixel_maximum_value=255;strcpy(name,"PhysBAM");colormap=0;}

    void Swap_Endian() const
    {Swap_Endianity(magic_number);Swap_Endianity(compression);Swap_Endianity(bytes_per_channel);Swap_Endianity(dimensions);
    Swap_Endianity(width);Swap_Endianity(height);Swap_Endianity(channels);Swap_Endianity(pixel_minimum_value);Swap_Endianity(pixel_maximum_value);Swap_Endianity(colormap);}        
};
}
#include <PhysBAM_Tools/Read_Write/Images/READ_WRITE_RGB_HEADER.h>
#endif
#endif
