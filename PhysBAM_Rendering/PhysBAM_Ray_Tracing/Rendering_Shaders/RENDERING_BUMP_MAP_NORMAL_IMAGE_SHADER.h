//#####################################################################
// Copyright 2004, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER  
//#####################################################################
#ifndef __RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER__
#define __RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_BUMP_MAP_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER:public RENDERING_BUMP_MAP_SHADER<T>
{
    typedef VECTOR<T,3> TV;
public:
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > normals;
    GRID<VECTOR<T,2> > grid;
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,VECTOR<T,3> > interpolation;
    
    RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER(const MATERIAL_SHADER<T>* shader,RENDER_WORLD<T>& world,const T s0_input=0,const T s1_input=1,const T t0_input=0,const T t1_input=1) 
        :RENDERING_BUMP_MAP_SHADER<T>(shader,world)
    {}

    virtual VECTOR<T,3> Perturbed_Normal(const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& object_space_position,const T s,const T t,const VECTOR<T,3>& object_space_normal,
        const VECTOR<T,3>& object_space_tangent,const VECTOR<T,3>& object_space_bitangent) const    
    {VECTOR<T,3> interpolated_normal=(interpolation.Clamped_To_Array_Cell(grid,normals,VECTOR<T,2>(s,t))).Normalized();
    MATRIX<T,3> tangent_space_to_object_space(object_space_tangent,object_space_bitangent,object_space_normal);
    return tangent_space_to_object_space*interpolated_normal;}

    void Initialize(const std::string& filename,const T max_phi=.01)
    {ARRAY<VECTOR<T,3> ,VECTOR<int,2> > pixels;IMAGE<T>::Read(filename,pixels);grid.Initialize(pixels.counts,RANGE<VECTOR<T,2> >::Unit_Box());
    normals.Resize(pixels.domain.min_corner.x,pixels.domain.max_corner.x,pixels.domain.min_corner.y,pixels.domain.max_corner.y);
    // convert each RGB value to its corresponding normal value
    for(int i=normals.domain.min_corner.x;i<=normals.domain.max_corner.x;++i)for(int j=normals.domain.min_corner.y;j<=normals.domain.max_corner.y;++j)normals(i,j)=pixels(i,j)*2-VECTOR<T,3>(1,1,1);}
    
    bool Valid() const
    {return grid.m>0 && grid.n>0 && grid.m==normals.m && grid.n==normals.n;}

//#####################################################################
};
}
#endif

