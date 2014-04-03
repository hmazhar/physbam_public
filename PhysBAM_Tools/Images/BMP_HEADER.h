//#####################################################################
// Copyright 2002-2005, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Struct BMP_HEADER
//#####################################################################
// A group of fields with the layout of a Windows bitmap file header. Used by the class BMP_FILE.
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __BMP_HEADER__
#define __BMP_HEADER__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
namespace PhysBAM{

struct BMP_HEADER{
    char file_type[2]; // must be "BM"
    int file_size;
    short reserved1;
    short reserved2;
    int offset;//#bytes from the beginning of the file to the bitmap
    int info_header_size;
    int w;
    int h;
    short number_of_bitplanes;//must be 1
    short bits_per_pixel;//e.g. 1,4,8, or 24
    int type_of_compression; //0 <=> no compression
    int bitmap_size;//w*h*3 if no compression
    int x_pixels_per_meter;
    int y_pixels_per_meter;
    int number_of_colors;
    int number_of_important_colors;

    void Initialize(const int w_input,const int h_input)
    {w=w_input;h=h_input;bitmap_size=w*h*3; //w*h*3 if no compression
    file_type[0]='B';file_type[1]='M';
    file_size=54+w*h*3;
    reserved1=0;reserved2=0;
    offset=54;//#bytes from the beginning of the file to the bitmap
    info_header_size=40;//probably unimportant
    number_of_bitplanes=1;//must be 1
    bits_per_pixel=24;//e.g. 1,4,8, or 24
    type_of_compression=0; //0 <=> no compression
    x_pixels_per_meter=0;//?
    y_pixels_per_meter=0;//?
    number_of_colors=0;
    number_of_important_colors=0;} 

//#####################################################################
};
}
#include <PhysBAM_Tools/Read_Write/Images/READ_WRITE_BMP_HEADER.h>
#endif
#endif
