//#####################################################################
// Copyright 2009-2010, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef COMPILE_WITHOUT_ZLIB_SUPPORT
#ifndef __ZIP__
#define __ZIP__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
#include <zlib.h>
namespace PhysBAM{
struct ZIP_FILE_HEADER;
//#####################################################################
// Functions Gzip_Out/Gzip_In - Create streams that read/write .gz
//#####################################################################
std::istream* Gzip_In(const std::string& filename,std::ios::openmode mode);
std::ostream* Gzip_Out(const std::string& filename,std::ios::openmode mode);
//#####################################################################
// Class ZIP_FILE_WRITER
//#####################################################################
class ZIP_FILE_WRITER
{
    std::ofstream ostream;
    ARRAY<ZIP_FILE_HEADER*> files;
public:

//#####################################################################
    ZIP_FILE_WRITER(const std::string& filename);
    virtual ~ZIP_FILE_WRITER();
    std::ostream* Add_File(const std::string& filename,const bool binary=true);
//#####################################################################
};

//#####################################################################
// Class ZIP_FILE_READER
//#####################################################################
class ZIP_FILE_READER
{
    std::ifstream istream;
public:
    HASHTABLE<std::string,ZIP_FILE_HEADER*> filename_to_header;
    
//#####################################################################
    ZIP_FILE_READER(const std::string& filename);
    virtual ~ZIP_FILE_READER();
    std::istream* Get_File(const std::string& filename,const bool binary=true);
    void Get_File_List(ARRAY<std::string>& filenames) const;
private:
    bool Find_And_Read_Central_Header();
//#####################################################################
};
}
#endif
#endif
#endif
