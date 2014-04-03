//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_ADAPTIVE
//#####################################################################
#ifndef __GRID_ADAPTIVE__
#define __GRID_ADAPTIVE__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
namespace PhysBAM{

template<class TV>
class GRID_ADAPTIVE:public GRID<TV>
{
    typedef GRID<TV> BASE;
    typedef typename TV::SCALAR T;    
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
public: 
    using BASE::Domain_Indices;
    //TODO: Use GRID_ADAPTIVE
    ARRAY<GRID<TV>*,TV_INT> sub_mac_grids;

    GRID_ADAPTIVE()
        :BASE()
    {Initialize_Sub_Grids();}

    GRID_ADAPTIVE(const int m_input,const int n_input,const int mn_input,const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const T zmin_input,const T zmax_input,
        const bool MAC_grid=false)
        :BASE(m_input,n_input,mn_input,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input,MAC_grid)
    {Initialize_Sub_Grids();}
    
    GRID_ADAPTIVE(const int m_input,const int n_input,const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const bool MAC_grid=false)
        :BASE(m_input,n_input,xmin_input,xmax_input,ymin_input,ymax_input,MAC_grid)
    {Initialize_Sub_Grids();}

    GRID_ADAPTIVE(const int m_input,const T xmin_input,const T xmax_input,const bool MAC_grid=false)
        :BASE(m_input,xmin_input,xmax_input,MAC_grid)
    {Initialize_Sub_Grids();}
    
    GRID_ADAPTIVE(const int m_input,const int n_input,const int mn_input,const RANGE<TV>& box,const bool MAC_grid=false)
        :BASE(m_input,n_input,mn_input,box,MAC_grid)
    {Initialize_Sub_Grids();}

    GRID_ADAPTIVE(const int m_input,const int n_input,const RANGE<TV>& box,const bool MAC_grid=false)
        :BASE(m_input,n_input,box,MAC_grid)
    {Initialize_Sub_Grids();}

    GRID_ADAPTIVE(const int m_input,const RANGE<TV>& box,const bool MAC_grid=false)
        :BASE(m_input,box,MAC_grid)
    {Initialize_Sub_Grids();}

    GRID_ADAPTIVE(const T dx,const RANGE<TV>& box,const bool MAC_grid=false)
        :BASE(dx,box,MAC_grid)
    {Initialize_Sub_Grids();}

    GRID_ADAPTIVE(const TV_INT& counts,const RANGE<TV>& box,const bool MAC_grid=false)
        :BASE(counts,box,MAC_grid)
    {Initialize_Sub_Grids();}

    GRID_ADAPTIVE(const GRID<TV>& grid_input)
        :BASE(grid_input)
    {Initialize_Sub_Grids();}

    void Initialize(const int m_input,const int n_input,const int mn_input,const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const T zmin_input,
        const T zmax_input,const bool MAC_grid=false)
    {Initialize(TV_INT(m_input,n_input,mn_input),RANGE<TV>(xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input),MAC_grid);}

    void Initialize(const int m_input,const int n_input,const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const bool MAC_grid=false)
    {Initialize(TV_INT(m_input,n_input),RANGE<TV>(xmin_input,xmax_input,ymin_input,ymax_input),MAC_grid);}

    void Initialize(const int m_input,const T xmin_input,const T xmax_input,const bool MAC_grid=false)
    {Initialize(TV_INT(m_input),RANGE<TV>(xmin_input,xmax_input),MAC_grid);}

    void Initialize(const int m_input,const int n_input,const int mn_input,const RANGE<TV>& box,const bool MAC_grid=false)
    {Initialize(TV_INT(m_input,n_input,mn_input),box,MAC_grid);}

    void Initialize(const int m_input,const int n_input,const RANGE<TV>& box,const bool MAC_grid=false)
    {Initialize(TV_INT(m_input,n_input),box,MAC_grid);}

    void Initialize(const int m_input,const RANGE<TV>& box,const bool MAC_grid=false)
    {Initialize(TV_INT(m_input),box,MAC_grid);}

    void Initialize(const T dx,const RANGE<TV>& box,const bool MAC_grid=false)
    {Initialize(TV_INT(box.Edge_Lengths()/dx),box,MAC_grid);}
    
    void Initialize(const TV_INT& counts_input,const RANGE<TV>& box,const bool MAC_grid=false)
    {BASE::Initialize(counts_input,box,MAC_grid);Initialize_Sub_Grids();}

    void Initialize_Sub_Grids()
    {sub_mac_grids.Resize(Domain_Indices());}
};
}
#endif


