//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header READ_WRITE_FORWARD
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_FORWARD__
#define __READ_WRITE_FORWARD__

namespace PhysBAM{

template<class T,class RW,class ENABLER=void> class Read_Write;
template<class T,class ENABLER=void> struct PLATFORM_INDEPENDENT_SIZE;
template<class T,class HAS=void> struct HAS_TYPED_READ;
template<class T,class HAS=void> struct HAS_TYPED_WRITE;

class STREAM_TYPE;
class TYPED_ISTREAM;
class TYPED_OSTREAM;

class PARAMETER_LIST;
template<class T> class IMAGE;
template<class T> class GENERIC_PARSER;

template<class T> class BMP_FILE;
template<class T> class EXR_FILE;
template<class T> class JPG_FILE;
template<class T> class PNG_FILE;
template<class T> class PPM_FILE;
template<class T> class RGB_FILE;

}
#endif
#endif
