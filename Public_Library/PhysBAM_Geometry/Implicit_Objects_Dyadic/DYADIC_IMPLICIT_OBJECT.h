//#####################################################################
// Copyright 2002-2007, Doug Enright, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Neil Molino, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DYADIC_IMPLICIT_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __DYADIC_IMPLICIT_SURFACE__
#define __DYADIC_IMPLICIT_SURFACE__

#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_POLICY.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_QUADTREE.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class TV>
class DYADIC_IMPLICIT_OBJECT:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename DYADIC_GRID_POLICY<TV>::DYADIC_GRID T_GRID;
    typedef typename INTERPOLATION_POLICY<T_GRID>::INTERPOLATION_SCALAR T_INTERPOLATION_SCALAR;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename T_GRID::CELL T_CELL;
    enum WORKAROUND {d=TV::m};
public:
    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;
    
    T_LEVELSET levelset;
    ARRAY<T>* phi_nodes;
    mutable ARRAY<VECTOR<T_CELL*,T_GRID::number_of_neighbors_per_cell> > neighbors;
private:
    T minimum_cell_size;
    bool need_destroy_data;
public:

    DYADIC_IMPLICIT_OBJECT(T_GRID& grid_input,ARRAY<T>& phi_input,ARRAY<T>* phi_nodes_input=0);
    ~DYADIC_IMPLICIT_OBJECT();

    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
    {levelset.grid.Fully_Refined_Block();}

    void Set_Custom_Interpolation(T_INTERPOLATION_SCALAR& interpolation_input)
    {levelset.Set_Custom_Interpolation(interpolation_input);}

    void Set_Custom_Secondary_Interpolation(T_INTERPOLATION_SCALAR& interpolation_input)
    {levelset.Set_Custom_Secondary_Interpolation(interpolation_input);}

    void Update_Box() PHYSBAM_OVERRIDE
    {box=levelset.grid.uniform_grid.domain;}

    void Update_Minimum_Cell_Size(const int maximum_depth=0) PHYSBAM_OVERRIDE
    {levelset.grid.Compute_Minimum_Cell_Constants(maximum_depth);
    minimum_cell_size=levelset.grid.Minimum_Edge_Length();}

    void Update_Phi_Nodes()
    {if(!phi_nodes) phi_nodes=new ARRAY<T>(levelset.grid.number_of_nodes,false);
    else phi_nodes->Resize(levelset.grid.number_of_nodes,false);
    LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::Interpolate_From_Cells_To_Nodes(levelset.grid,levelset.phi,*phi_nodes);}
    
    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const PHYSBAM_OVERRIDE
    {T result=levelset.grid.uniform_grid.Minimum_Edge_Length();
    ARRAY<T_CELL*> cells_intersecting_box;levelset.grid.Get_Cells_Intersecting_Box(box,cells_intersecting_box);
    for(int i=1;i<=cells_intersecting_box.m;i++)result=min(result,cells_intersecting_box(i)->DX().Min());return result;}

    T operator()(const TV& location) const PHYSBAM_OVERRIDE
    {return levelset.Phi(location,phi_nodes);}

    T operator()(const BLOCK_DYADIC<T_GRID>& block,const TV& location) const
    {return levelset.Phi(block,location);}

    T Extended_Phi(const TV& location) const PHYSBAM_OVERRIDE
    {return levelset.Extended_Phi(location,phi_nodes);}

    T Phi_Secondary(const TV& location) const PHYSBAM_OVERRIDE
    {return levelset.Phi_Secondary(location,phi_nodes);}

    TV Normal(const TV& location,const int aggregate=-1) const PHYSBAM_OVERRIDE
    {if(aggregate != -1) return box.Normal(aggregate);else return levelset.Normal(location,true,phi_nodes);}

    TV Extended_Normal(const TV& location,const int aggregate=-1) const PHYSBAM_OVERRIDE
    {if(aggregate != -1) return box.Normal(aggregate);else return levelset.Extended_Normal(location,true,phi_nodes);}

    SYMMETRIC_MATRIX<T,d> Hessian(const TV& X) const PHYSBAM_OVERRIDE
    {return levelset.Hessian(X);}

    VECTOR<T,d-1> Principal_Curvatures(const TV& X) const PHYSBAM_OVERRIDE
    {return levelset.Principal_Curvatures(X);}
    
    void Rescale(const T scaling_factor) PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}// This function will most definitely not work without a significant effort (and cost)

    void Inflate(const T inflation_distance) PHYSBAM_OVERRIDE
    {levelset.phi-=inflation_distance;}

    T Integration_Step(const T phi) const PHYSBAM_OVERRIDE
    {return (T).1*max(abs(phi),minimum_cell_size);}

    T Minimum_Cell_Size() const  PHYSBAM_OVERRIDE
    {return minimum_cell_size;}

    T Min_Phi() const  PHYSBAM_OVERRIDE
    {return levelset.phi.Min();}

    bool Lazy_Inside(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return box.Lazy_Inside(location) && levelset.Lazy_Inside(location,contour_value,phi_nodes);}
    
    bool Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return box.Lazy_Inside(location) && levelset.Lazy_Inside_And_Value(location,phi_value,contour_value,phi_nodes);}

    bool Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return levelset.Lazy_Inside_Extended_Levelset(unclamped_X,contour_value,phi_nodes);}

    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return levelset.Lazy_Inside_Extended_Levelset_And_Value(unclamped_X,phi_value,contour_value,phi_nodes);}
    
    bool Lazy_Outside(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return box.Lazy_Outside(location) || levelset.Lazy_Outside(location,contour_value,phi_nodes);}

    bool Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return levelset.Lazy_Outside_Extended_Levelset(unclamped_X,contour_value,phi_nodes);}

    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return levelset.Lazy_Outside_Extended_Levelset_And_Value(unclamped_X,phi_value,contour_value,phi_nodes);}

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name()
    {return STRING_UTILITIES::string_sprintf("DYADIC_IMPLICIT_OBJECT<T,VECTOR<T,%d> >",TV::dimension);}

    virtual std::string Extension() const PHYSBAM_OVERRIDE {return Static_Extension();}
    static std::string Static_Extension()
    {return TV::dimension==2?"quad":"oct";}

    bool Intersection(RAY<TV>& ray,const T thickness=0) const PHYSBAM_OVERRIDE
    {return Intersection(ray,thickness,*phi_nodes);}
    
//###########################################################################
    static DYADIC_IMPLICIT_OBJECT* Create();
    bool Intersection(RAY<TV>& ray,const T thickness,const ARRAY<T>& phi_nodes,const bool verbose=false) const;
//###########################################################################
};
}
#endif
#endif
