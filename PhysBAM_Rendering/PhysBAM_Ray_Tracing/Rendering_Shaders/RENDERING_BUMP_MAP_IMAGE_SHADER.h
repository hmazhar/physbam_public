//#####################################################################
// Copyright 2004, Geoffrey Irving, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_BUMP_MAP_IMAGE_SHADER  
//#####################################################################
#ifndef __RENDERING_BUMP_MAP_IMAGE_SHADER__
#define __RENDERING_BUMP_MAP_IMAGE_SHADER__

#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_BUMP_MAP_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_BUMP_MAP_IMAGE_SHADER:public RENDERING_BUMP_MAP_SHADER<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<T,2> TV2;
public:
    GRID<TV2> grid;
    ARRAY<T,VECTOR<int,2> > phi;
    LEVELSET_2D<GRID<TV2> > levelset;
    T s0,s_min,s_max,s_scaling_factor,t0,t_min,t_max,t_scaling_factor;
    
    RENDERING_BUMP_MAP_IMAGE_SHADER(const MATERIAL_SHADER<T>* shader,RENDER_WORLD<T>& world,const T s0_input=0,const T s1_input=1,const T t0_input=0,const T t1_input=1) 
        :RENDERING_BUMP_MAP_SHADER<T>(shader,world),levelset(grid,phi),s0(s0_input),s_min(min(s0_input,s1_input)),s_max(max(s0_input,s1_input)),s_scaling_factor((T)1/(s1_input-s0_input)),
        t0(t0_input),t_min(min(t0_input,t1_input)),t_max(max(t0_input,t1_input)),t_scaling_factor((T)1/(t1_input-t0_input))
    {}
    
    VECTOR<T,3> Perturbed_Normal(const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& object_space_position,const T s,const T t,const VECTOR<T,3>& object_space_normal,const VECTOR<T,3>& object_space_tangent,const VECTOR<T,3>& object_space_bitangent) const PHYSBAM_OVERRIDE
    {if(s<s_min || s>s_max || t<t_min || t>t_max)return object_space_normal;VECTOR<T,2> location((s-s0)*s_scaling_factor,(t-t0)*t_scaling_factor);
    T fx=(levelset.Phi(VECTOR<T,2>(location.x+grid.dX.x,location.y))-levelset.Phi(VECTOR<T,2>(location.x-grid.dX.x,location.y)))/(2*grid.dX.x),
        fy=(levelset.Phi(VECTOR<T,2>(location.x,location.y+grid.dX.y))-levelset.Phi(VECTOR<T,2>(location.x,location.y-grid.dX.y)))/(2*grid.dX.y);
    VECTOR<T,3> N0(0,0,1),N1(-fx,-fy,1); return (object_space_normal+(N1-N0)).Normalized();}
    
    bool Initialize(const std::string& filename,const T max_phi=.01)
    {ARRAY<VECTOR<T,3> ,VECTOR<int,2> > pixels;IMAGE<T>::Read(filename,pixels);
    int i,j;//for(i=pixels.domain.min_corner.x;i<=pixels.domain.max_corner.x;i++)for(j=pixels.domain.min_corner.y;j<pixels.domain.min_corner.y+pixels.domain.max_corner.y-j;j++) exchange(pixels(i,j),pixels(i,pixels.domain.max_corner.y+pixels.domain.min_corner.y-j));
    grid.Initialize(pixels.counts,RANGE<VECTOR<T,2> >::Unit_Box()); 
    phi.Resize(pixels.domain);
    for(i=phi.domain.min_corner.x;i<=phi.domain.max_corner.x;++i)for(j=phi.domain.min_corner.y;j<=phi.domain.max_corner.y;++j)phi(i,j)=VECTOR<T,3>::Dot_Product(VECTOR<T,3>((T).299,(T).587,(T).114),pixels(i,j));
    phi*=max_phi/phi.Max();
    return true;}
    
    bool Valid() const
    {return grid.m>0 && grid.n>0 && grid.m==phi.m && grid.n==phi.n;}

    void Print() const
    {LOG::cout<<"grid "<<grid<<std::endl;}

//#####################################################################
};
}
#endif

