//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_LIGHT__
#define __RENDERING_LIGHT__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
namespace PhysBAM{

template<class T> class RENDER_WORLD;
template<class T> class RENDERING_RAY;
template<class T> class PHOTON_MAP;

template<class T>
class RENDERING_LIGHT:public NONCOPYABLE
{
public:
    int light_index;
    VECTOR<T,3> position;
    VECTOR<T,3> color;
    T brightness;
    RENDER_WORLD<T>& world;
    bool supports_global_photon_mapping,supports_caustic_photon_mapping,supports_volume_photon_mapping,photon_source_only,casts_shadows;
    mutable RANDOM_NUMBERS<T> sample_points_random,global_photon_random,caustic_photon_random,volume_photon_random;

    RENDERING_LIGHT(const VECTOR<T,3>& position_input,const VECTOR<T,3>& color_input,const T brightness_input,RENDER_WORLD<T>& world_input,const bool supports_global_photons=false,
        const bool supports_caustic_photons=false,const bool supports_volume_photons=false,const bool photon_source=false)
        :light_index(-1),position(position_input),color(color_input),brightness(brightness_input),world(world_input),supports_global_photon_mapping(supports_global_photons),
        supports_caustic_photon_mapping(supports_caustic_photons),supports_volume_photon_mapping(supports_volume_photons),photon_source_only(photon_source),casts_shadows(true)
    {
        global_photon_random.Set_Seed(1);caustic_photon_random.Set_Seed(2);volume_photon_random.Set_Seed(3);sample_points_random.Set_Seed(4);
    }

    virtual ~RENDERING_LIGHT()
    {}

//#####################################################################
    virtual int Emit_Photons(RENDERING_RAY<T>& parent_ray,PHOTON_MAP<T>& photon_map,const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type) const {return 0;}
    virtual void Sample_Points(const VECTOR<T,3>& surface_position,const VECTOR<T,3>& surface_normal,ARRAY<RAY<VECTOR<T,3> > >& sample_array) const=0;
    virtual VECTOR<T,3> Emitted_Light(const RENDERING_RAY<T>& ray) const=0;
//#####################################################################
};
}
#endif
