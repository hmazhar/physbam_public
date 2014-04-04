//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Duc Nguyen, Andrew Selle, Justin Solomon.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RAY_TRACING_DRIVER  
//##################################################################### 
#ifndef __RAY_TRACING_DRIVER__
#define __RAY_TRACING_DRIVER__

#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Tools/Math_Tools/integer_log.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include "RAY_TRACING_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class RAY_TRACING_DRIVER
{
public:
    RAY_TRACING_EXAMPLE<T>& example;
    RENDER_WORLD<T> world;
    PROGRESS_INDICATOR progress; 
    std::string basedir;
    int large_pixel_size;
    std::string output_filename,alpha_filename;

    RAY_TRACING_DRIVER(RAY_TRACING_EXAMPLE<T>& example_input)
        :example(example_input)
    {}

    virtual ~RAY_TRACING_DRIVER()
    {}

//#####################################################################
    void Initialize();
    virtual void Execute_Main_Program();
//#####################################################################
};
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void RAY_TRACING_DRIVER<T>::
Initialize()
{
    // Call back to example to setup
    LOG::Time("Initialize Scene");
    example.Initialize_Scene(world,example.frame);

    // determine pixels to render
    RANGE<VECTOR<int,2> >& clip=example.clipping_region;
    clip.min_corner.x=max(clip.min_corner.x,1);clip.max_corner.x=min(clip.max_corner.x,world.camera.film.grid.counts.x);
    clip.min_corner.y=max(clip.min_corner.y,1);clip.max_corner.y=min(clip.max_corner.y,world.camera.film.grid.counts.y);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    large_pixel_size=world.use_irradiance_cache?1<<integer_log(clip.max_corner.x-clip.min_corner.x):1;
#endif

    output_filename=example.Get_Output_Filename(example.frame);
    alpha_filename=example.Get_Alpha_Filename(example.frame);
    basedir=FILE_UTILITIES::Get_Base_Directory_Name(output_filename);FILE_UTILITIES::Create_Directory(basedir); // initialize output directory
    if(!IMAGE<T>::Is_Supported(output_filename)){LOG::cerr<<"Image format for '"<<output_filename<<"' not supported!"<<std::endl;exit(1);}
    if(FILE_UTILITIES::File_Exists(output_filename)){if(!FILE_UTILITIES::File_Writable(output_filename)){LOG::cerr<<"Output file '"<<output_filename<<"' not writable!"<<std::endl;exit(1);}}
    else if(!FILE_UTILITIES::Directory_Writable(basedir)){LOG::cerr<<"Directory '"<<basedir<<"' not writable!"<<std::endl;exit(1);}

    // TODO: make sure alpha file is writable

    LOG::Time("Prepare for Forward Ray Trace");
    world.Prepare_For_Forward_Ray_Trace();
}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
template<class T> void RAY_TRACING_DRIVER<T>::
Execute_Main_Program()
{    
    LOG::SCOPE scope("RENDER","Rendering frame %d",example.frame);
    Initialize();
    // render
    LOG::Time("Forward Ray Trace");
    world.Render(example.clipping_region, progress);

    //TODO: Integrate something similar to this code back in so that we can ray trace with a preview while using best-candidate sampling
    //BOX<VECTOR<int,2> >& clip=example.clipping_region;
    //ARRAY<VECTOR<T,2> > sample_positions;ARRAY<VECTOR<T,3> > sample_values;
    //CAMERA<T>& camera=world.camera;
    //FILM<T>& film=camera.film;
    //ARRAY<T,VECTOR<int,2> > sum_weights(1,film.grid.m,1,film.grid.n);ARRAY<VECTOR<T,3> ,VECTOR<int,2> > rgb_ghost(film.grid,3);
    //for(int pixel_size=large_pixel_size;pixel_size>0;pixel_size/=2) for(int i=clip.min_corner.x;i<=clip.max_corner.x;i++)for(int j=clip.min_corner.y;j<=clip.max_corner.y;j++)
    //    if(i%pixel_size==0 && j%pixel_size==0 && (pixel_size==large_pixel_size || i%(pixel_size*2)!=0 || j%(pixel_size*2)!=0)){
    //        // put the ray color on film
    //        //world.Set_Pixel_Color(i,j,world.Compute_Pixel_Color(i,j));
    //        sample_positions.Remove_All();sample_values.Remove_All();
    //        world.Render_Pixel(VECTOR<int,2>(i,j),0);          
    //        // print a percentage progress if requested
    //        progress.Progress(); 
    //        if(example.pixels_between_output_directory_existence_checks && progress.done%example.pixels_between_output_directory_existence_checks==0
    //            && !FILE_UTILITIES::Directory_Exists(basedir)){
    //            LOG::cout<<"RAY_TRACING_DRIVER<T>::Execute_Main_Program(): base directory "<<basedir<<" doesn't exist!"<<std::endl;exit(1);}}

    // write images
    LOG::Time("Writing Images");
    world.camera.film.Print_Film_Clipped(output_filename,example.gamma_correction,example.clipping_region);
}
//#####################################################################
}
#endif
