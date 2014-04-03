//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// This is mainly from the book "Advanced Renderman"
//#####################################################################
#ifndef __RENDERING_WOOD_SHADER__
#define __RENDERING_WOOD_SHADER__

#include <PhysBAM_Tools/Random_Numbers/NOISE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_WOOD_SHADER:public MATERIAL_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:

    VECTOR<T,3> color1,color2;

    T ring_frequency;
    T ring_noise,ring_noise_frequency;
    T trunk_wobble,trunk_wobble_frequency;
    T angular_wobble,angular_wobble_frequency;
    T grain_frequency;
    T grainy;
    T ringy;

    T Smooth(const T value,const T min,const T max) const
    {T v=clamp((value-min)/(max-min),(T)0,(T)1);
    return v*v*(3-2*v);}

    T Smooth_Pulse(T e0,T e1,T e2,T e3,T x) const
    {return Smooth(x,e0,e1)-Smooth(x,e2,e3);}

    T Smooth_Pulse_Train(T e0,T e1,T e2,T e3,T period,T x) const
    {double ipart;
    T mod=(T)modf(x/period,&ipart);
    mod*=period;
    return Smooth_Pulse(e0,e1,e2,e3,mod);}

    RENDERING_WOOD_SHADER(const VECTOR<T,3>& color1_input,const VECTOR<T,3>& color2_input,const T ring_frequency,const T ring_noise,const T ring_noise_frequency,const T trunk_wobble,
        const T trunk_wobble_frequency,const T angular_wobble,const T angular_wobble_frequency,const T grain_frequency,const T grainy,const T ringy,RENDER_WORLD<T>& world) 
        :MATERIAL_SHADER<T>(world),color1(color1_input),color2(color2_input),ring_frequency(ring_frequency),ring_noise(ring_noise),ring_noise_frequency(ring_noise_frequency),
        trunk_wobble(trunk_wobble),trunk_wobble_frequency(trunk_wobble_frequency),angular_wobble(angular_wobble),angular_wobble_frequency(angular_wobble_frequency),
        grain_frequency(grain_frequency),grainy(grainy),ringy(ringy)
    {}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {RAY<VECTOR<T,3> > object_space_ray=intersection_object.Object_Space_Ray(ray.ray);VECTOR<T,3> object_space_position=object_space_ray.Point(object_space_ray.t_max);
    TV p=intersection_object.Get_Solid_Texture_Coordinates(object_space_position,ray.ray.aggregate_id);
    TV Pshad=p;
    // wobble the rings
    TV offset=NOISE<T>::Fractional_Brownian_3D(Pshad*ring_noise_frequency,2,4,(T).5);
    TV Pring=Pshad+ring_noise*offset;
    // make the trunk not quite on the z-axis
    Pring+=trunk_wobble*NOISE<T>::Noise_Function3(TV(0,0,Pshad.z*trunk_wobble_frequency))*TV(1,1,0);
    // evaluate the radius
    T r=Pring.Remove_Index(3).Magnitude()*ring_frequency;
    // add some noise to the angle to simulate rings that are not quite round
    r+=angular_wobble*Smooth(r,(T)0,(T)5)*NOISE<T>::Noise_Function1(angular_wobble_frequency*Pring*TV(1,1,(T).1));
    r+=(T)0.5*NOISE<T>::Noise_Function1(TV(r,0,0));
    T in_ring=Smooth_Pulse_Train((T).1,(T).55,(T).7,(T).95,(T)1,r);
    // now the grain
    TV Pgrain=Pshad*grain_frequency*TV(1,1,(T).05);
    T grain=0;
    T amp=1;
    T dPgrain=(T).01;
    for(int i=0;i<2;i++){
        T grainvalid=1-Smooth((T).4,(T).2,dPgrain);
        if(grainvalid){
            T g=grainvalid*NOISE<T>::Noise_Function1(Pgrain);
            g*=((T).3+(T).7*in_ring);
            g=pow(clamp((T).8-g,(T)0,(T)1),(T)2);
            g=grainy*Smooth(g,(T)0.5,(T)1);
            if(i==0) in_ring*=(1-(T).4*grainvalid);
            grain=max(grain,g);}
        dPgrain*=2;
        Pgrain*=2;
        amp*=(T).5;}
    // blend between colors
    T a=(1-grain)*in_ring*ringy+grain;
    return (1-a)*color1+a*color2;}

//#####################################################################
};
}
#endif
