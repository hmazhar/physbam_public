//#####################################################################
// Copyright 2005, Jiayi Chong, Jeong-Mo Hong, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_LEVELSET_MULTIPLE_OBJECT
//#####################################################################
#ifndef __RENDERING_LEVELSET_MULTIPLE_OBJECT__
#define __RENDERING_LEVELSET_MULTIPLE_OBJECT__

#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_LEVELSET_MULTIPLE_REGION_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>

namespace PhysBAM{

template<class T_LEVELSET_MULTIPLE>
class RENDERING_LEVELSET_MULTIPLE_OBJECT:public RENDERING_OBJECT<typename T_LEVELSET_MULTIPLE::VECTOR_T::SCALAR>
{
    typedef typename T_LEVELSET_MULTIPLE::VECTOR_T TV;typedef typename TV::SCALAR T;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::inverse_transform;using RENDERING_OBJECT<T>::priority;

    T_LEVELSET_MULTIPLE levelset_multiple;
    ARRAY<RENDERING_LEVELSET_MULTIPLE_REGION_OBJECT<T,T_LEVELSET_MULTIPLE>* > rendering_levelset_multiple_region_objects;
    int number_of_regions;
 
    RENDERING_LEVELSET_MULTIPLE_OBJECT(GRID<TV>& grid_input,ARRAY<ARRAY<T,VECTOR<int,3> > >& phis_input)        
        :levelset_multiple(grid_input,phis_input),number_of_regions(phis_input.m)
    {rendering_levelset_multiple_region_objects.Resize(phis_input.m);
    for(int i=1;i<=number_of_regions;i++) rendering_levelset_multiple_region_objects(i)=new RENDERING_LEVELSET_MULTIPLE_REGION_OBJECT<T,T_LEVELSET_MULTIPLE>(levelset_multiple,i);}

    virtual ~RENDERING_LEVELSET_MULTIPLE_OBJECT()
    {for(int i=1;i<rendering_levelset_multiple_region_objects.m;i++) delete rendering_levelset_multiple_region_objects(i);}

    bool Intersection(RAY<VECTOR<T,3> >& ray,const int lowest_priority,const RENDERING_OBJECT<T>** intersected_object) const PHYSBAM_OVERRIDE
    {if(priority<lowest_priority) return false;int region_start,region_end;        
    if(Intersection(ray,region_start,region_end,small_number)){
        if(region_end==-1) *intersected_object=rendering_levelset_multiple_region_objects(region_start);
        else if(region_start==-1) *intersected_object=rendering_levelset_multiple_region_objects(region_end);
        else{
            if(region_start>0&&rendering_levelset_multiple_region_objects(region_start)->priority>rendering_levelset_multiple_region_objects(region_end)->priority)
                *intersected_object=rendering_levelset_multiple_region_objects(region_start);
            else *intersected_object=rendering_levelset_multiple_region_objects(region_end);}
        return true;}
    return false;}

    int Intersected_Region(RAY<VECTOR<T,3> >& ray) const
    {int region_start,region_end;        
    if(Intersection(ray,region_start,region_end,small_number)){
        if(region_end==-1) return region_start;
        else if(region_start==-1) return region_end;
        else{
            if(region_start>0&&rendering_levelset_multiple_region_objects(region_start)->priority>rendering_levelset_multiple_region_objects(region_end)->priority)
                return region_start;
            else return region_end;}
        return -1;}
    return -1;}

    bool Inside_Region_Only(const VECTOR<T,3>& location,int region_check) const
    {bool is_inside=false;
      for(int i=1;i<=number_of_regions;i++) {
        bool inside_region=rendering_levelset_multiple_region_objects(i)->Inside(location);
        if(i==region_check && inside_region) is_inside=true;
        if(i!=region_check && inside_region) return false;
      }
      return true;}

    bool Inside(const VECTOR<T,3>& location,RENDERING_OBJECT<T>** intersected_object) const PHYSBAM_OVERRIDE
    {for(int i=1;i<=number_of_regions;i++) if(rendering_levelset_multiple_region_objects(i)->Inside(location)){
        *intersected_object=(RENDERING_OBJECT<T>*)rendering_levelset_multiple_region_objects(i);return true;}
    return false;}

    RANGE<VECTOR<T,3> > Object_Space_Bounding_Box() const  PHYSBAM_OVERRIDE
    {return levelset_multiple.grid.domain;}

    T Integration_Step(const T phi) const
    {T distance=abs(phi);
    if(distance > 3*levelset_multiple.grid.min_dX) return (T).5*distance;    
    else if(distance > levelset_multiple.grid.min_dX) return (T).25*distance;
    return (T).1*levelset_multiple.grid.min_dX;}

    TRIANGULATED_SURFACE<T>* Generate_Triangles()const PHYSBAM_OVERRIDE {return TRIANGULATED_SURFACE<T>::Create();};

//#####################################################################
    bool Intersection(RAY<VECTOR<T,3> >& ray,int& region_start,int& region_end,const T thickness=0) const;
//#####################################################################
};   
}
#endif

