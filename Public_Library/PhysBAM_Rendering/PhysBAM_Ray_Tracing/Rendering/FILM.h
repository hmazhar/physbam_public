//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Joyce Pan, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FILM
//#####################################################################
#ifndef __FILM__
#define __FILM__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Functions.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{
template<class T> class CAMERA;

template<class T> T Film_Box_Filter(const VECTOR<T,2>& X,const VECTOR<T,2>& width)
{if(abs(X.x)<width.x && abs(X.y)<width.y) return 1;else return 0;}

template<class T> T Film_Gaussian_Filter(const VECTOR<T,2>& X,const VECTOR<T,2>& width)
{return exp(-2*(sqr(2*X.x/width.x)+sqr(2*X.y/width.y)));}

template<class T> T Lanczos_1D(T x,T tau)
{if(abs(x)>(T)1) return (T)0;x*=(T)pi;return sinc(x*tau)*sinc(x);}

template<class T> T Film_Lanczos_Filter(const VECTOR<T,2>& X,const VECTOR<T,2>& width)
{T tau=(T).5;return Lanczos_1D(X.x/width.x,tau)*Lanczos_1D(X.y/width.y,tau);}

template<class T>
class FILM
{
    typedef VECTOR<T,2> TV2;typedef VECTOR<T,3> TV;
    
public:
    struct SAMPLE{
        T time;
        VECTOR<T,2> film_position;
        VECTOR<T,3> world_position;
        VECTOR<T,3> radiance;
        T alpha; // for now, probably either 0 or 1
    };
    typedef T (*PIXEL_FILTER)(const VECTOR<T,2>& pixel_X,const VECTOR<T,2>& filter_width);

    T width,height; // size of the film - input from the camera
    GRID<TV2> grid; // 2D grid of pixel locations
    GRID<TV2> subpixel_grid; // 2D grid of pixel locations
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > colors; // color at each pixel
    ARRAY<T,VECTOR<int,2> > weights; // filter weights
    ARRAY<T,VECTOR<int,2> > alphas; // percent of pixel occupied
    int samples_per_pixel;
    bool use_four_subpixels;
    RANDOM_NUMBERS<T> random;
    VECTOR<T,2> filter_width;
    VECTOR<T,2> effective_filter_width;
    PIXEL_FILTER filter;
    T dither_amplitude;
    static const int best_candidate_sample_count;
    static const T best_candidate_samples[4096][2];

    enum SAMPLER{BEST_CANDIDATE,UNIFORM,JITTERED};
    SAMPLER sampler;

    FILM();
    ~FILM();

    void Set_Filter(const VECTOR<T,2>& width_input,PIXEL_FILTER filter_input)
    {filter=filter_input;filter_width=width_input;effective_filter_width=grid.dX*filter_width;}

    void Print_Film(const std::string& filename,const T gamma) const
    {Print_Film_Clipped(filename,gamma,grid.Domain_Indices());}

    TV Pixel_Color(const VECTOR<int,2>& pixel) const
    {return weights(pixel.x,pixel.y)*colors(pixel.x,pixel.y);}

    VECTOR<int,2> Samples_Extent()
    {if(sampler==UNIFORM || sampler==JITTERED) return VECTOR<int,2>(1,1);
        else /*if sampler==BEST_CANDIDATE*/ return (int)sqrt((T)best_candidate_sample_count/(T)samples_per_pixel)*VECTOR<int,2>(1,1);}

//#####################################################################
    void Set_Resolution(const int pixels_width,const int pixels_height);
    void Print_Film_Clipped(const std::string& filename,const T gamma,const RANGE<VECTOR<int,2> >& box) const;
    void Generate_Samples(const RANGE<VECTOR<int,2> >& pixel_range,const CAMERA<T>& camera,ARRAY<SAMPLE>& samples);
    void Generate_Stratified_Sample_Points(const VECTOR<int,2>& pixel_index,const CAMERA<T>& camera,ARRAY<SAMPLE>& samples);
    void Add_Sample(const SAMPLE& sample);
//#####################################################################
};
}

#endif
