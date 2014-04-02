//#####################################################################
// Copyright 2003-2006, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEEP_WATER_EVOLUTION
//#####################################################################
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/DEEP_WATER_EVOLUTION.h>
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/FREQUENCY_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void DEEP_WATER_EVOLUTION<TV>::
Initialize()
{
    random.Set_Seed(12345);

    if(grid.counts/2*2!=grid.counts){
        LOG::cerr<<"Deep water requires even numbers of grid cells (grid = "<<grid<<")"<<std::endl;
        PHYSBAM_FATAL_ERROR();}
    h.Resize(grid.Domain_Indices());
    if(lambda){ // displacements for choppy waves
        for(int k=1;k<=TV::m;k++) displacement[k].Resize(grid.Domain_Indices());
        Xh.Resize(grid.Domain_Indices());}

    h_hat.Resize(FREQUENCY_ITERATOR<TV>::Domain_Indices(grid));
    h_hat1.Resize(FREQUENCY_ITERATOR<TV>::Domain_Indices(grid));
    h_hat2.Resize(FREQUENCY_ITERATOR<TV>::Domain_Indices(grid));
    for(int axis=1;axis<=TV::m;axis++) dXh_hat[axis].Resize(FREQUENCY_ITERATOR<TV>::Domain_Indices(grid));

    // Initialize FFT stuff
    fft.grid=grid;

    if(phillips_spectrum.amplitude){
        T scale=sqr((T)grid.counts.Product()); // scale Fourier coefficients to cancel the scaling that will occur in Inverse_Transform (need sqr since amplitude isn't a normal amplitude)
        phillips_spectrum.grid=grid;
        phillips_spectrum.gravity=gravity;
        phillips_spectrum.amplitude*=scale;
        phillips_spectrum.Generate_Spectrum(random);
        phillips_spectrum.amplitude/=scale;
        phillips_spectrum.Get_Separate_H_Hats(h_hat1,h_hat2);
        phillips_spectrum.h_initial.Clean_Memory();}

    if(TV::dimension==2 && texture_cutoffs!=TV_INT())
        for(FREQUENCY_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT I=iterator.Index();
            if((I.x<texture_cutoffs.x || I.x>(grid.counts.x-texture_cutoffs.x)) && I[2]<texture_cutoffs[2])
                h_hat1(I)=h_hat2(I)=COMPLEX<T>();}

    Advance_Height(0);
}
//#####################################################################
// Function Set_H_Hats_From_Height
//#####################################################################
template<class TV> void DEEP_WATER_EVOLUTION<TV>::
Set_H_Hats_From_Height()
{
    fft.Transform(h,h_hat);
    for(FREQUENCY_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT I=iterator.Index();
        COMPLEX<T> delta=h_hat(I)-h_hat1(I)-h_hat2(I);
        h_hat1(I)+=(T).5*delta;
        h_hat2(I)+=(T).5*delta;}

    // set DC to 0
    h_hat1(TV_INT())=h_hat2(TV_INT())=COMPLEX<T>();
}
//#####################################################################
// Function Advance_H_Hats
//#####################################################################
template<class TV> void DEEP_WATER_EVOLUTION<TV>::
Advance_H_Hats(const T dt)
{
    T_ARRAYS_COMPLEX scaled_p_hat;
    if(use_surface_pressure){
        scaled_p_hat.Resize(FREQUENCY_ITERATOR<TV>::Domain_Indices(grid));
        fft.Transform(surface_pressure,scaled_p_hat);scaled_p_hat/=gravity;
        for(FREQUENCY_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT I=iterator.Index();
            COMPLEX<T> half_delta=(T).5*(scaled_p_hat(I)-h_hat1(I)-h_hat2(I));
            h_hat1(I)+=half_delta;h_hat2(I)+=half_delta;}}

    // dispersion relation
    for(FREQUENCY_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT I=iterator.Index();
        T k=iterator.Frequency().Magnitude(),omega=Omega(k);
        COMPLEX<T> exp_dt=COMPLEX<T>::Unit_Polar(omega*dt);

        h_hat1(I)*=exp_dt;
        h_hat2(I)*=exp_dt.Conjugated();}
}
//#####################################################################
// Function Advance_Height
//#####################################################################
template<class TV> void DEEP_WATER_EVOLUTION<TV>::
Advance_Height(const T dt)
{
    Advance_H_Hats(dt);

    // Initialize the (0) elements since they won't get set below
    for(int axis=1;axis<=TV::m;axis++) dXh_hat[axis](TV_INT())=COMPLEX<T>();

    for(FREQUENCY_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT I=iterator.Index();
        TV k1=iterator.Frequency();
        T k=k1.Magnitude();
        h_hat(I)=h_hat1(I)+h_hat2(I);
        if(k) for(int axis=1;axis<=TV::m;axis++)
            dXh_hat[axis](I)=COMPLEX<T>(0,-k1[axis]/k)*h_hat(I);}

    // Set DC to 0
    h_hat(TV_INT())=COMPLEX<T>();

    if(filter_high_frequencies) fft.Filter_High_Frequencies(h_hat,high_frequency_cutoff);
    fft.Inverse_Transform(h_hat,h,true,false);
    if(lambda){
        for(int k=1;k<=TV::m;k++) fft.Inverse_Transform(dXh_hat[k],displacement[k],true,false);
        for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            TV DX;for(int axis=1;axis<=TV::m;axis++) DX[axis]=displacement[axis](node);
            Xh(node)=grid.X(node)+lambda*DX;}}
}
//#####################################################################
// Function Get_Vertical_Velocity
//#####################################################################
template<class TV> void DEEP_WATER_EVOLUTION<TV>::
Get_Vertical_Velocity(T_ARRAYS_T& v) const
{
    T_ARRAYS_COMPLEX v_hat(FREQUENCY_ITERATOR<TV>::Domain_Indices(grid));
    for(FREQUENCY_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT I=iterator.Index();
        T k=iterator.Frequency().Magnitude(),omega=Omega(k);
        v_hat(I)=COMPLEX<T>(0,omega)*(h_hat1(I)-h_hat2(I));}
    fft.Inverse_Transform(v_hat,v,true,false);
}
//#####################################################################
// Function Texture_Shallow_Water
//#####################################################################
template<class TV> void DEEP_WATER_EVOLUTION<TV>::
Texture_Shallow_Water(const int frame,const STREAM_TYPE stream_type,const std::string& shallow_water_directory)
{
    std::string f=FILE_UTILITIES::Number_To_String(frame);
    T_ARRAYS_T shallow_water_eta;FILE_UTILITIES::Read_From_File(stream_type,shallow_water_directory+"/eta."+f,shallow_water_eta);
    T_ARRAYS_T shallow_water_ground;FILE_UTILITIES::Read_From_File(stream_type,shallow_water_directory+"/ground",shallow_water_eta);
    for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
        T shallow_water_h=shallow_water_eta(node)-shallow_water_ground(node);
        T min_h=1;
        T max_h=10;
        T alpha=clamp((shallow_water_h-min_h)/(max_h-min_h),(T)0,(T)1);
        h(node)=shallow_water_eta(node)+alpha*h(node);
#if 0
        if(shallow_water_eta(node)-shallow_water_ground(node)<1e-5)
            h(node)=shallow_water_eta(node);
        else{
            if(h(node)+shallow_water_eta(node)<shallow_water_ground(node))
                h(node)=shallow_water_eta(node);
            else
                h(node)+=shallow_water_eta(node);}
#endif
    }
}
//#####################################################################
template class DEEP_WATER_EVOLUTION<VECTOR<float,1> >;
template class DEEP_WATER_EVOLUTION<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEEP_WATER_EVOLUTION<VECTOR<double,1> >;
template class DEEP_WATER_EVOLUTION<VECTOR<double,2> >;
#endif
