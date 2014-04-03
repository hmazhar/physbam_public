//#####################################################################
// Copyright 2008, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Images/MOV_FILE.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#ifdef USE_LIBJPEG
extern "C"{
#include <jpeglib.h>
}
#endif
namespace PhysBAM{

typedef unsigned int uint;
typedef unsigned short ushort;

namespace{
    //#####################################################################
    // Function Big_Endian
    //#####################################################################
    template<class T>
    static void Big_Endian(T& x)
    {
        assert(sizeof(T)<=8);
        if(sizeof(T)>1) {T old=x;for(unsigned int k=1;k<=sizeof(T);k++) ((char*)&x)[k-1]=((char*)&old)[sizeof(T)-k];}
    }
    //#####################################################################
    // Function Write
    //#####################################################################
    static void Write(FILE* fp, char num)
    {
        fwrite(&num,sizeof(num),1,fp);
    }
    //#####################################################################
    // Function Write
    //#####################################################################
    static void Write(FILE* fp, ushort num)
    {
#ifndef PHYSBAM_BIG_ENDIAN
        Big_Endian(num);
#endif
        fwrite(&num,sizeof(num),1,fp);
    }
    //#####################################################################
    // Function Wrte
    //#####################################################################
    static void Write(FILE* fp, uint num)
    {
#ifndef PHYSBAM_BIG_ENDIAN
        Big_Endian(num);
#endif
        fwrite(&num,sizeof(num),1,fp);
    }
    //#####################################################################
    // Function Write_Identity_Matrix
    //#####################################################################
    static void Write_Identity_Matrix(FILE* fp)
    {
        Write(fp,(uint)0x10000);Write(fp,(uint)0x00000);Write(fp,(uint)0); // 16.16 fixed pt
        Write(fp,(uint)0x00000);Write(fp,(uint)0x10000);Write(fp,(uint)0); // 16.16 fixed pt
        Write(fp,(uint)0x00000);Write(fp,(uint)0x00000);Write(fp,(uint)0x40000000); // 2.30 fixed pt
    }
}

class QT_ATOM
{
    FILE *fp;
    long start_offset;
    const char* type;
public:
    QT_ATOM(FILE* fp,const char* type)
        :fp(fp),type(type)
    {
        start_offset=ftell(fp);
        uint dummy;
        fwrite(&dummy,4,1,fp);
        fputs(type,fp);
    }

    ~QT_ATOM()
    {
        uint atom_size=ftell(fp)-start_offset;
        Big_Endian(atom_size);
        fseek(fp,start_offset,SEEK_SET);
        fwrite(&atom_size,4,1,fp);
        fseek(fp,0,SEEK_END);
    }

    long Offset()
    {return start_offset;}
    //#####################################################################
};


//#####################################################################
// Function MOV_WRTER
//#####################################################################
template<class T> MOV_WRITER<T>::
MOV_WRITER(const std::string& filename,const int frames_per_second)
    :frames_per_second(frames_per_second),width(0),height(0)
{
    fp=fopen(filename.c_str(),"wb");
    if(!fp) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Failed to open %s for writing",filename.c_str()));
    current_mov=new QT_ATOM(fp,"mdat");
}
//#####################################################################
// Function ~MOV_WRTER
//#####################################################################
template<class T> MOV_WRITER<T>::
~MOV_WRITER()
{
    delete current_mov;
    Write_Footer();
    fclose(fp);
}
//#####################################################################
// Function Add_Frame
//#####################################################################
template<class T> void MOV_WRITER<T>::
Add_Frame(ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image)
{
#ifdef USE_LIBJPEG
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    if(width==0 && height==0){width=image.counts.x;height=image.counts.y;}
    if(width!=image.counts.x || height!=image.counts.y) throw std::runtime_error("Frame does not have same size as previous frame(s)");
    
    cinfo.err=jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    long frame_begin=ftell(fp);
    jpeg_stdio_dest(&cinfo,fp     );
    cinfo.image_width=image.counts.x;cinfo.image_height=image.counts.y;cinfo.input_components=3;
    cinfo.in_color_space=JCS_RGB; // colorspace of input image
    jpeg_set_defaults(&cinfo);jpeg_set_quality(&cinfo,95,TRUE); // limit to baseline-JPEG values
    jpeg_start_compress(&cinfo,TRUE);
    
    int row_stride=cinfo.image_width*3; // JSAMPLEs per row in image_buffer
    JSAMPLE* row=new unsigned char[row_stride];JSAMPROW row_pointer[]={row};
    while(cinfo.next_scanline < cinfo.image_height){
        int index=0;for(int i=1;i<=image.counts.x;i++){VECTOR<unsigned char,3> pixel=IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,image.counts.y-cinfo.next_scanline));row[index++]=pixel.x;row[index++]=pixel.y;row[index++]=pixel.z;} // copy row
        jpeg_write_scanlines(&cinfo,row_pointer,1);}
    delete[] row;
    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);
    long frame_end=ftell(fp);
    sample_lengths.Append(frame_end-frame_begin);
    sample_offsets.Append(frame_begin-current_mov->Offset());
#endif
}
//#####################################################################
// Function Write_Footer
//#####################################################################
template<class T> void MOV_WRITER<T>::
Write_Footer()
{
    const int frames=sample_offsets.m;assert(sample_offsets.m==sample_lengths.m);
    QT_ATOM a(fp,"moov");
        {QT_ATOM a(fp,"mvhd");
            char c=0;
            Write(fp,c); // version
            Write(fp,c);Write(fp,c);Write(fp,c); // reserved
            Write(fp,(uint)0); // creation time
            Write(fp,(uint)0); // modification time
            Write(fp,(uint)frames_per_second); // time rate 1/30th second
            Write(fp,(uint)frames); // duration
            Write(fp,(uint)0x10000); // preferred rate 16bit fixed pt 1.0
            Write(fp,(ushort)0x100); // full volume
            Write(fp,(uint)0);Write(fp,(uint)0);Write(fp,(ushort)0x0); // 10 bytes padded
            Write_Identity_Matrix(fp);
            Write(fp,(uint)0); // preview time
            Write(fp,(uint)0); // preview duration
            Write(fp,(uint)0); // poster time 
            Write(fp,(uint)0); // selection time
            Write(fp,(uint)frames); // selection duration
            Write(fp,(uint)0); // current time
            Write(fp,(uint)2);} // next track
        {QT_ATOM a(fp,"trak");
            {QT_ATOM a(fp,"tkhd");
                Write(fp,(uint)0xf); // flag visibble
                Write(fp,(uint)0); // creation time
                Write(fp,(uint)0); // modification time
                Write(fp,(uint)1); // track id
                Write(fp,(uint)0); // reserved
                Write(fp,(uint)frames); // duration
                Write(fp,(uint)0);Write(fp,(uint)0); // reserved
                Write(fp,(ushort)0); // layer
                Write(fp,(ushort)0); // alternative group
                Write(fp,(ushort)0x100); // volume
                Write(fp,(ushort)0); // reserved
                Write_Identity_Matrix(fp);
                Write(fp,(uint)width<<16); // width
                Write(fp,(uint)height<<16);} // height
            {QT_ATOM a(fp,"edts");
              {QT_ATOM a(fp,"elst");
                Write(fp,(uint)0); // version flags
                Write(fp,(uint)0);}} // 1 entry
            {QT_ATOM a(fp,"mdia");
                {QT_ATOM a(fp,"mdhd");
                    Write(fp,(uint)0x0); // version/flag visibble
                    Write(fp,(uint)0); // creation time
                    Write(fp,(uint)0); // modified time
                    Write(fp,(uint)frames_per_second); // time scale
                    Write(fp,(uint)frames); // duration
                    Write(fp,(ushort)0); // english language
                    Write(fp,(ushort)0xffff);} // quality
                {QT_ATOM a(fp,"hdlr");
                    Write(fp,(uint)0x0); // version/flags
                    fputs("mhlrvide",fp);
                    Write(fp,(uint)0); // component manufacture
                    Write(fp,(uint)0); // component flags
                    Write(fp,(uint)0); // component flags mask
                    Write(fp,(char)0); // component name
                    fputs("Linux Video Media Handler",fp);}
                {QT_ATOM a(fp,"minf");
                    {QT_ATOM a(fp,"vmhd");
                        Write(fp,(uint)0x0001); // version/flags set 1 for compatibility
                        Write(fp,(ushort)0x40); // graphics mode copy
                        Write(fp,(ushort)0x8000); // unused graphics mode opcolor
                        Write(fp,(ushort)0x8000); // unused graphics mode opcolor
                        Write(fp,(ushort)0x8000);} // unused graphics mode opcolor
                    {QT_ATOM a(fp,"hdlr");
                        Write(fp,(uint)0x0); // version/flags
                        fputs("dhlralis",fp);
                        Write(fp,(uint)0); // component manufacture
                        Write(fp,(uint)0); // component flags
                        Write(fp,(uint)0); // component flags mask
                        Write(fp,(char)0); // component name
                        fputs("Linux Alias Data Handler",fp);}
                    {QT_ATOM a(fp,"dinf");
                        {QT_ATOM a(fp,"dref");
                            Write(fp,(uint)0x0); // vvvf version flags
                            Write(fp,(uint)0x1); // 1 entry
                            {QT_ATOM a(fp,"alis");
                                Write(fp,(uint)1);}}}
                    {QT_ATOM a(fp,"stbl");
                        {QT_ATOM a(fp,"stsd");
                            Write(fp,(uint)0); // version and flags
                            Write(fp,(uint)1); // 1 entry
                            {QT_ATOM a(fp,"jpeg");
                                Write(fp,(uint)0); //reserved
                                Write(fp,(ushort)0); //reserved
                                Write(fp,(ushort)1); // data reference index
                                // write video specific data
                                Write(fp,(ushort)0); // version
                                Write(fp,(ushort)0); // revision level
                                fputs("lnux",fp); // vendor
                                Write(fp,(uint)100); //temporal quality (max)
                                Write(fp,(uint)258); //spatial quality (max)
                                Write(fp,(ushort)width); // width of image
                                Write(fp,(ushort)height); // height of image
                                Write(fp,(uint)0x00480000); // height of image (72dpi)
                                Write(fp,(uint)0x00480000); // height of image (72dpi)
                                Write(fp,(uint)0); // data size (must be zero)
                                Write(fp,(ushort)1); // frames per sample (usually 1)
                                const char* descript="Quicktime for Linux";
                                Write(fp,(char)strlen(descript)); 
                                assert(strlen(descript)<32);
                                fputs(descript,fp); // compressor
                                for(size_t i=0;i<32-strlen(descript)-1;i++) Write(fp,(char)0);
                                Write(fp,(ushort)24); // color depth
                                Write(fp,(ushort)(short)-1);}} // use default color table id
                        {QT_ATOM a(fp,"stts");
                            Write(fp,(uint)0); // version and flags
                            Write(fp,(uint)1); // 1 entry
                            Write(fp,(uint)frames); // all frames have same duration
                            Write(fp,(uint)1);} // duration is one time unit
                        {QT_ATOM a(fp,"stsc");
                            Write(fp,(uint)0); // version and flags
                            Write(fp,(uint)1); // 1 entry
                            Write(fp,(uint)1); // first sample to use in chunk
                            Write(fp,(uint)1); // number of samples per chunk
                            Write(fp,(uint)1);} // index of descriptor (points to stsd)
                        {QT_ATOM a(fp,"stsz");
                            Write(fp,(uint)0); // version and flags
                            Write(fp,(uint)0); // sample size (non-uniform so zero and table follows)
                            Write(fp,(uint)frames); // one entry per frame
                            for(int i=1;i<=sample_lengths.Size();i++) Write(fp,(uint)sample_lengths(i));}
                        {QT_ATOM a(fp,"stco");
                            Write(fp,(uint)0); // version and flags
                            Write(fp,(uint)frames); // one entry per frame
                            for(int i=1;i<=sample_offsets.Size();i++) Write(fp,(uint)sample_offsets(i));}}}}} // offset from begin of file
}
//#####################################################################
// Function Enabled
//#####################################################################
template<class T> bool MOV_WRITER<T>::
Enabled()
{
#ifdef USE_LIBJPEG
    return true; 
#else
    return false;
#endif
}
//#####################################################################
template class MOV_WRITER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MOV_WRITER<double>;
#endif
}
#endif
