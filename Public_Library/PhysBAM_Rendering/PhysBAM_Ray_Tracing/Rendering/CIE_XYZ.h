//#####################################################################
// Copyright 2002, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CIE_XYZ
//#####################################################################
#ifndef __CIE_XYZ__
#define __CIE_XYZ__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class CIE_XYZ
{
    typedef VECTOR<T,1> TV;
public:
    int tablesize;
    T starting_wavelength,ending_wavelength; // in meters
    GRID<TV> grid; // in meters
    ARRAY<T,VECTOR<int,1> > X_spectrum,Y_spectrum,Z_spectrum; // from the CIE tables
    T display_adaption_luminance;  // default is 50 candelas/m^2
    T display_maximum_luminance; // default is 100 candelas/m^2
    T world_adaption_luminance; // needs to be calculated and set in units of candelas/m^2
    T world_to_display; // display_luminance=world_to_display*world_luminance
    T world_to_unit;      // unit_luminance=world_to_unit*world_luminance;
    
    CIE_XYZ()
        :tablesize(81),starting_wavelength((T)380e-9),ending_wavelength((T)780e-9),
        grid(tablesize,starting_wavelength,ending_wavelength),
        X_spectrum(1,grid.counts.x),Y_spectrum(1,grid.counts.x),Z_spectrum(1,grid.counts.x)
    {      
        Initialize_XYZ_From_CIE_Tables();
        Set_Display_Adaption_Luminance();
        Set_Display_Maximum_Luminance();
    }

    void Set_Display_Adaption_Luminance(const T display_adaption_luminance_input=50) // units are candela/m^2
    {display_adaption_luminance=display_adaption_luminance_input;}

    void Set_Display_Maximum_Luminance(const T display_maximum_luminance_input=100) // units are candela/m^2
    {display_maximum_luminance=display_maximum_luminance_input;}
    
    void Set_World_Adaption_Luminance(const T world_adaption_luminance_input) // units are candela/m^2
    {world_adaption_luminance=world_adaption_luminance_input;}

//#####################################################################
    void Initialize_XYZ_From_CIE_Tables();
    T Radiometric_To_Photometric(const ARRAY<T,VECTOR<int,1> >& radiometric_spectrum) const; 
    VECTOR<T,3> Calculate_XYZ(const ARRAY<T,VECTOR<int,1> >& radiometric_radiance_spectrum) const; 
    void Calculate_Luminance_Scale_Factors();
    VECTOR<T,3> World_To_Display_XYZ(const VECTOR<T,3>& XYZ) const;
    VECTOR<T,3> World_To_Unit_XYZ(const VECTOR<T,3>& XYZ) const;
    VECTOR<T,3> XYZ_To_RGB(const VECTOR<T,3>& XYZ) const;
    VECTOR<T,3> RGB_To_XYZ(const VECTOR<T,3>& RGB) const;
//#####################################################################
};   
}
#endif
