//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Igor Neverov, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace FILE_UTILITIES
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __FILE_UTILITIES__
#define __FILE_UTILITIES__

#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <PhysBAM_Tools/Utilities/EXCEPTIONS.h>
namespace PhysBAM{

namespace FILE_UTILITIES{

//#####################################################################
// ADD NEW FILE EXTENSIONS HERE
// The enumeration should match the file_extensions array, and
// UNKNOWN_FILE should be matched with a 0 in the array.
enum FILE_TYPE{RGD_FILE,RGD2D_FILE,TRI_FILE,PHI_FILE,PHI2D_FILE,OCT_FILE,PLY_FILE,PLY2D_FILE,RGB_FILE,TRI2D_FILE,CURVE_FILE,CURVE2D_FILE,TET_FILE,HEX_FILE,BOX_FILE,PHONEME_FILE,UNKNOWN_FILE};
//#####################################################################

//###################################################################
// Platform Specific Function Definitions
//###################################################################
bool Directory_Exists(const std::string& dirname);
bool Create_Directory(const std::string& dirname);
std::string Real_Path(const std::string& path);
int Compare_File_Times_Ignoring_Compression_Suffix(const std::string& filename1,const std::string& filename2);
std::istream* Safe_Open_Input(const std::string& filename,bool binary=true);
std::ostream* Safe_Open_Output(const std::string& filename,bool binary=true,bool write_compressed_if_possible=true);
FILE* Temporary_File();
//###################################################################
// Platform Non-Specific Function Declarations
//###################################################################
int Compare_File_Times(const std::string& filename1,const std::string& filename2);
bool File_Exists_Ignoring_Compression_Suffix(const std::string& filename);
bool File_Writable_Ignoring_Compression_Suffix(const std::string& filename);
bool Directory_Writable(const std::string& dirname);
std::string Get_Working_Directory();
//###################################################################
// Platform Non-Specific Function Definitions
//###################################################################
std::string Get_File_Extension_Ignoring_Compression_Suffix(const std::string &filename);
std::string Get_Basename_Ignoring_Compression_Suffix(const std::string& filename);
std::string Get_Short_Name_Ignoring_Compression_Suffix(const std::string& filename);
bool File_Extension_Matches_Ignoring_Compression_Suffix(const std::string& filename,const std::string& ext,const bool case_sensitive=true);
bool File_Is_Compressed(const std::string& filename);
bool File_Exists(const std::string& filename);
bool File_Writable(const std::string& filename);
std::string Strip_Compression_Suffix(const std::string& filename);
std::string Real_File(const std::string& filename);
bool File_Extension_Matches(const std::string& filename,const std::string& ext,const bool case_sensitive=true);
FILE_TYPE Get_File_Type_Ignoring_Compression_Suffix(const std::string& filename);
FILE_TYPE Get_File_Type(const std::string& filename);
bool File_Type_Matches_Ignoring_Compression_Suffix(const std::string& filename,FILE_TYPE type);
bool File_Type_Matches(const std::string& filename,FILE_TYPE type);

inline bool Is_Rgd_File(const std::string& filename)
{return File_Type_Matches(filename,RGD_FILE);}

inline bool Is_Rgd2D_File(const std::string& filename)
{return File_Type_Matches(filename,RGD2D_FILE);}

inline bool Is_Tri_File(const std::string& filename)
{return File_Type_Matches(filename,TRI_FILE);}

inline bool Is_Tet_File(const std::string& filename)
{return File_Type_Matches(filename,TET_FILE);}

inline bool Is_Hex_File(const std::string& filename)
{return File_Type_Matches(filename,HEX_FILE);}

inline bool Is_Phi_File(const std::string& filename)
{return File_Type_Matches(filename,PHI_FILE);}

inline bool Is_Phi2D_File(const std::string& filename)
{return File_Type_Matches(filename,PHI2D_FILE);}

inline bool Is_Oct_File(const std::string& filename)
{return File_Type_Matches(filename,OCT_FILE);}

inline bool Is_Ply_File(const std::string& filename)
{return File_Type_Matches(filename,PLY_FILE);}

inline bool Is_Ply2D_File(const std::string& filename)
{return File_Type_Matches(filename,PLY2D_FILE);}

inline bool Is_Rgb_File(const std::string& filename)
{return File_Type_Matches(filename,RGB_FILE);}

inline bool Is_Curve_File(const std::string& filename)
{return File_Type_Matches(filename,CURVE_FILE);}

inline bool Is_Curve2D_File(const std::string& filename)
{return File_Type_Matches(filename,CURVE2D_FILE);}

inline bool Is_Box_File(const std::string& filename)
{return File_Type_Matches(filename,BOX_FILE);}

inline bool Is_Phoneme_File(const std::string& filename)
{return File_Type_Matches(filename,PHONEME_FILE);}

std::string Get_File_Extension(const std::string &filename);
std::string Get_Basename(const std::string& filename);
std::string Get_Base_Directory_Name(const std::string& path);
std::string Get_Short_Name(const std::string& filename);

inline std::string Number_To_String(const int i)
{return STRING_UTILITIES::Value_To_String(i);}

std::string Find_First_Nonexistent_File_In_Sequence(std::string filename_pattern,const int id_start=0,int* id_result=0);
std::string Find_First_Nonexistent_Directory_In_Sequence(std::string directory_pattern,const int id_start=0,int* id_final=0);
std::string Make_First_Nonexistent_Directory_In_Sequence(std::string directory_pattern,const int id_start=0,int* id_final=0);
//###################################################################
// Utilities for "animated" files (files with %d for frame number)
//###################################################################
bool Is_Animated(const std::string &filename);
std::string Get_Frame_Filename(const std::string &filename,int frame);
bool Frame_File_Exists(const std::string &filename,int frame);
//#####################################################################
// Read_From_File
//#####################################################################
// Convenience functions
template<class RW,class T1>
inline void Read_From_File(const std::string& filename,T1& d1)
{std::istream* input=Safe_Open_Input(filename);Read_Binary<RW>(*input,d1);delete input;}
template<class RW,class T1,class T2>
inline void Read_From_File(const std::string& filename,T1& d1,T2& d2)
{std::istream* input=Safe_Open_Input(filename);Read_Binary<RW>(*input,d1,d2);delete input;}
template<class RW,class T1,class T2,class T3>
inline void Read_From_File(const std::string& filename,T1& d1,T2& d2,T3& d3)
{std::istream* input=Safe_Open_Input(filename);Read_Binary<RW>(*input,d1,d2,d3);delete input;}
template<class RW,class T1,class T2,class T3,class T4>
inline void Read_From_File(const std::string& filename,T1& d1,T2& d2,T3& d3,T4& d4)
{std::istream* input=Safe_Open_Input(filename);Read_Binary<RW>(*input,d1,d2,d3,d4);delete input;}
template<class RW,class T1,class T2,class T3,class T4,class T5>
inline void Read_From_File(const std::string& filename,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5)
{std::istream* input=Safe_Open_Input(filename);Read_Binary<RW>(*input,d1,d2,d3,d4,d5);delete input;}
template<class RW,class T1,class T2,class T3,class T4,class T5,class T6>
inline void Read_From_File(const std::string& filename,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6)
{std::istream* input=Safe_Open_Input(filename);Read_Binary<RW>(*input,d1,d2,d3,d4,d5,d6);delete input;}
template<class RW,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline void Read_From_File(const std::string& filename,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6,T7& d7)
{std::istream* input=Safe_Open_Input(filename);Read_Binary<RW>(*input,d1,d2,d3,d4,d5,d6,d7);delete input;}
template<class RW,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
inline void Read_From_File(const std::string& filename,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6,T7& d7,T8& d8)
{std::istream* input=Safe_Open_Input(filename);Read_Binary<RW>(*input,d1,d2,d3,d4,d5,d6,d7,d8);delete input;}
// runtime float/double versions
template<class T1>
inline void Read_From_File(const STREAM_TYPE stream_type,const std::string& filename,T1& d1)
{std::istream* input=Safe_Open_Input(filename);TYPED_ISTREAM typed_input(*input,stream_type);Read_Binary(typed_input,d1);delete input;}
template<class T1,class T2>
inline void Read_From_File(const STREAM_TYPE stream_type,const std::string& filename,T1& d1,T2& d2)
{std::istream* input=Safe_Open_Input(filename);TYPED_ISTREAM typed_input(*input,stream_type);Read_Binary(typed_input,d1,d2);delete input;}
template<class T1,class T2,class T3>
inline void Read_From_File(const STREAM_TYPE stream_type,const std::string& filename,T1& d1,T2& d2,T3& d3)
{std::istream* input=Safe_Open_Input(filename);TYPED_ISTREAM typed_input(*input,stream_type);Read_Binary(typed_input,d1,d2,d3);delete input;}
template<class T1,class T2,class T3,class T4>
inline void Read_From_File(const STREAM_TYPE stream_type,const std::string& filename,T1& d1,T2& d2,T3& d3,T4& d4)
{std::istream* input=Safe_Open_Input(filename);TYPED_ISTREAM typed_input(*input,stream_type);Read_Binary(typed_input,d1,d2,d3,d4);delete input;}
template<class T1,class T2,class T3,class T4,class T5>
inline void Read_From_File(const STREAM_TYPE stream_type,const std::string& filename,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5)
{std::istream* input=Safe_Open_Input(filename);TYPED_ISTREAM typed_input(*input,stream_type);Read_Binary(typed_input,d1,d2,d3,d4,d5);delete input;}
template<class T1,class T2,class T3,class T4,class T5,class T6>
inline void Read_From_File(const STREAM_TYPE stream_type,const std::string& filename,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6)
{std::istream* input=Safe_Open_Input(filename);TYPED_ISTREAM typed_input(*input,stream_type);Read_Binary(typed_input,d1,d2,d3,d4,d5,d6);delete input;}
template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline void Read_From_File(const STREAM_TYPE stream_type,const std::string& filename,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6,T7& d7)
{std::istream* input=Safe_Open_Input(filename);TYPED_ISTREAM typed_input(*input,stream_type);Read_Binary(typed_input,d1,d2,d3,d4,d5,d6,d7);delete input;}
template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
inline void Read_From_File(const STREAM_TYPE stream_type,const std::string& filename,T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6,T7& d7,T8& d8)
{std::istream* input=Safe_Open_Input(filename);TYPED_ISTREAM typed_input(*input,stream_type);Read_Binary(typed_input,d1,d2,d3,d4,d5,d6,d7,d8);delete input;}
//#####################################################################
// Write_To_File
//#####################################################################
// Convenience functions
template<class RW,class T1>
inline void Write_To_File(const std::string& filename,const T1& d1)
{std::ostream* output=Safe_Open_Output(filename);Write_Binary<RW>(*output,d1);delete output;}
template<class RW,class T1,class T2>
inline void Write_To_File(const std::string& filename,const T1& d1,const T2& d2)
{std::ostream* output=Safe_Open_Output(filename);Write_Binary<RW>(*output,d1,d2);delete output;}
template<class RW,class T1,class T2,class T3>
inline void Write_To_File(const std::string& filename,const T1& d1,const T2& d2,const T3& d3)
{std::ostream* output=Safe_Open_Output(filename);Write_Binary<RW>(*output,d1,d2,d3);delete output;}
template<class RW,class T1,class T2,class T3,class T4>
inline void Write_To_File(const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4)
{std::ostream* output=Safe_Open_Output(filename);Write_Binary<RW>(*output,d1,d2,d3,d4);delete output;}
template<class RW,class T1,class T2,class T3,class T4,class T5>
inline void Write_To_File(const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5)
{std::ostream* output=Safe_Open_Output(filename);Write_Binary<RW>(*output,d1,d2,d3,d4,d5);delete output;}
template<class RW,class T1,class T2,class T3,class T4,class T5,class T6>
inline void Write_To_File(const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6)
{std::ostream* output=Safe_Open_Output(filename);Write_Binary<RW>(*output,d1,d2,d3,d4,d5,d6);delete output;}
template<class RW,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline void Write_To_File(const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6,const T7& d7)
{std::ostream* output=Safe_Open_Output(filename);Write_Binary<RW>(*output,d1,d2,d3,d4,d5,d6,d7);delete output;}
template<class RW,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
inline void Write_To_File(const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6,const T7& d7,const T8& d8)
{std::ostream* output=Safe_Open_Output(filename);Write_Binary<RW>(*output,d1,d2,d3,d4,d5,d6,d7,d8);delete output;}
// runtime float/double versions
template<class T1>
inline void Write_To_File(const STREAM_TYPE stream_type,const std::string& filename,const T1& d1)
{std::ostream* output=Safe_Open_Output(filename);TYPED_OSTREAM typed_output(*output,stream_type);Write_Binary(typed_output,d1);delete output;}
template<class T1,class T2>
inline void Write_To_File(const STREAM_TYPE stream_type,const std::string& filename,const T1& d1,const T2& d2)
{std::ostream* output=Safe_Open_Output(filename);TYPED_OSTREAM typed_output(*output,stream_type);Write_Binary(typed_output,d1,d2);delete output;}
template<class T1,class T2,class T3>
inline void Write_To_File(const STREAM_TYPE stream_type,const std::string& filename,const T1& d1,const T2& d2,const T3& d3)
{std::ostream* output=Safe_Open_Output(filename);TYPED_OSTREAM typed_output(*output,stream_type);Write_Binary(typed_output,d1,d2,d3);delete output;}
template<class T1,class T2,class T3,class T4>
inline void Write_To_File(const STREAM_TYPE stream_type,const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4)
{std::ostream* output=Safe_Open_Output(filename);TYPED_OSTREAM typed_output(*output,stream_type);Write_Binary(typed_output,d1,d2,d3,d4);delete output;}
template<class T1,class T2,class T3,class T4,class T5>
inline void Write_To_File(const STREAM_TYPE stream_type,const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5)
{std::ostream* output=Safe_Open_Output(filename);TYPED_OSTREAM typed_output(*output,stream_type);Write_Binary(typed_output,d1,d2,d3,d4,d5);delete output;}
template<class T1,class T2,class T3,class T4,class T5,class T6>
inline void Write_To_File(const STREAM_TYPE stream_type,const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6)
{std::ostream* output=Safe_Open_Output(filename);TYPED_OSTREAM typed_output(*output,stream_type);Write_Binary(typed_output,d1,d2,d3,d4,d5,d6);delete output;}
template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline void Write_To_File(const STREAM_TYPE stream_type,const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6,const T7& d7)
{std::ostream* output=Safe_Open_Output(filename);TYPED_OSTREAM typed_output(*output,stream_type);Write_Binary(typed_output,d1,d2,d3,d4,d5,d6,d7);delete output;}
template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
inline void Write_To_File(const STREAM_TYPE stream_type,const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6,const T7& d7,const T8& d8)
{std::ostream* output=Safe_Open_Output(filename);TYPED_OSTREAM typed_output(*output,stream_type);Write_Binary(typed_output,d1,d2,d3,d4,d5,d6,d7,d8);delete output;}
//#####################################################################
// Read_From_Text_File
//#####################################################################
// Convenience function
template<class T1>
inline void Read_From_Text_File(const std::string& filename,T1& d1)
{std::istream* input=Safe_Open_Input(filename,false);*input>>d1;delete input;}
//#####################################################################
// Write_To_Text_File
//#####################################################################
// Convenience functions
template<class T1>
inline void Write_To_Text_File(const std::string& filename,const T1& d1)
{std::ostream* output=Safe_Open_Output(filename,false);*output<<d1;delete output;}
template<class T1,class T2>
inline void Write_To_Text_File(const std::string& filename,const T1& d1,const T2& d2)
{std::ostream* output=Safe_Open_Output(filename,false);*output<<d1<<d2;delete output;}
template<class T1,class T2,class T3>
inline void Write_To_Text_File(const std::string& filename,const T1& d1,const T2& d2,const T3& d3)
{std::ostream* output=Safe_Open_Output(filename,false);*output<<d1<<d2<<d3;delete output;}
template<class T1,class T2,class T3,class T4>
inline void Write_To_Text_File(const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4)
{std::ostream* output=Safe_Open_Output(filename,false);*output<<d1<<d2<<d3<<d4;delete output;}
template<class T1,class T2,class T3,class T4,class T5>
inline void Write_To_Text_File(const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5)
{std::ostream* output=Safe_Open_Output(filename,false);*output<<d1<<d2<<d3<<d4<<d5;delete output;}
template<class T1,class T2,class T3,class T4,class T5,class T6>
inline void Write_To_Text_File(const std::string& filename,const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6)
{std::ostream* output=Safe_Open_Output(filename,false);*output<<d1<<d2<<d3<<d4<<d5<<d6;delete output;}
//#####################################################################
// Create_From_File
//#####################################################################
template<class T>
inline void Create_From_File(const STREAM_TYPE stream_type,const std::string& filename,T*& d)
{typename REMOVE_CONST<T>::TYPE* read=T::Create();Read_From_File(stream_type,filename,*read);d=read;}

template<class RW,class T>
inline void Create_From_File(const std::string& filename,T*& d)
{Create_From_File(STREAM_TYPE(RW()),filename,d);}

template<class T>
inline T* Create_From_File(const STREAM_TYPE stream_type,const std::string& filename)
{typename REMOVE_CONST<T>::TYPE* d=T::Create();Read_From_File(stream_type,filename,*d);return d;}

void Ignore(std::istream& input,char c);
//#####################################################################
}
}
#endif
#endif
