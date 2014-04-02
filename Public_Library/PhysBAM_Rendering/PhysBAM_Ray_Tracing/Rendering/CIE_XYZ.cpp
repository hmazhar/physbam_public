//#####################################################################
// Copyright 2002, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//##################################################################### 
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/CIE_XYZ.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize_XYZ_From_CIE_Tables
//#####################################################################
// from the CIE tables at http://www.hike.te.chiba-u.ac.jp/ikeda/CIE/home.html
template<class T> void CIE_XYZ<T>::
Initialize_XYZ_From_CIE_Tables() 
{
    int i;
    
    // CIE 15.2-1986 Table 2.1 data
    // CIE 1931 Standard Colorimetric Observer x2(lambda) data between 380 nm and 780 nm at 5 nm intervals
    double X_table[]={.0014,.0022,.0042,.0077,.0143,.0232,.0435,.0776,.1344,.2148,.2839,.3285,.3483,.3481,.3362,.3187,.2908,
                                .2511,.1954,.1421,.0956,.0580,.0320,.0147,.0049,.0024,.0093,.0291,.0633,.1096,.1655,.2257,.2904,.3597,
                                .4334,.5121,.5945,.6784,.7621,.8425,.9163,.9786,1.0263,1.0567,1.0622,1.0456,1.0026,.9384,.8544,.7514,
                                .6424,.5419,.4479,.3608,.2835,.2187,.1649,.1212,.0874,.0636,.0468,.0329,.0227,.0158,.0114,.0081,.0058,
                                .0041,.0029,.0020,.0014,.0010,.0007,.0005,.0003,.0002,.0002,.0001,.0001,.0001,0};
    for(i=1;i<=81;i++) X_spectrum(i)=(T)X_table[i-1];

    // CIE 15.2-1986 Table 2.1 data
    // CIE 1931 Standard Colorimetric Observer y2(lambda) data between 380 nm and 780 nm at 5 nm intervals
    double Y_table[]={0,.0001,.0001,.0002,.0004,.0006,.0012,.0022,.004,.0073,.0116,.0168,.023,.0298,.038,.048,.06,.0739,.091,
                                .1126,.139,.1693,.208,.2586,.323,.4073,.503,.6082,.71,.7932,.862,.9149,.954,.9803,.995,1,.995,.9786,.952,
                                .9154,.87,.8163,.757,.6949,.631,.5668,.503,.4412,.381,.321,.265,.217,.175,.1382,.107,.0816,.061,.0446,
                                .032,.0232,.017,.0119,.0082,.0057,.0041,.0029,.0021,.0015,.001,.0007,.0005,.0004,.0002,.0002,.0001,
                                .0001,.0001,0,0,0,0};
    for(i=1;i<=81;i++) Y_spectrum(i)=(T)Y_table[i-1];
    
    // CIE 15.2-1986 Table 2.1 data
    // CIE 1931 Standard Colorimetric Observer z2(lambda) data between 380 nm and 630 nm at 5 nm intervals
    double Z_table[]={.0065,.0105,.0201,.0362,.0679,.1102,.2074,.3713,.6456,1.0391,1.3856,1.623,1.7471,1.7826,1.7721,
                                1.7441,1.6692,1.5281,1.2876,1.0419,.813,.6162,.4652,.3533,.272,.2123,.1582,.1117,.0782,.0573,.0422,
                                .0298,.0203,.0134,.0087,.0057,.0039,.0027,.0021,.0018,.0017,.0014,.0011,.001,.0008,.0006,.0003,.0002,
                                .0002,.0001,0};
    for(i=1;i<=51;i++) Z_spectrum(i)=(T)Z_table[i-1];
}
//#####################################################################
// Function Radiometric_To_Photometric
//#####################################################################
// assumes the radiometric_spectrum is defined on the grid from 380 nm to 780 nm in 5 nm intervals
// energy_spectrum - from joules/m to talbots
// power_spectrum - from watts/m to lumens
// intensity_spectrum - from from watts/steradian/m to candela
// radiosity_specturm - from watts/m^2/m to lumens/m^2
// radiance_spectrum - from watts/(steradian*m^2)/m to candela/m^2
template<class T> T CIE_XYZ<T>::
Radiometric_To_Photometric(const ARRAY<T,VECTOR<int,1> >& radiometric_spectrum) const
{
    T conversion_constant=680; // units are lumens/watt
    T photometric_scalar=0;
    for(int i=1;i<=grid.counts.x;i++) photometric_scalar+=Y_spectrum(i)*radiometric_spectrum(i);
    photometric_scalar*=(grid.dX.x*conversion_constant);
    return photometric_scalar;
}
//#####################################################################
// Function Calculate_XYZ
//#####################################################################
// assumes the radiometric_radiance is defined on the grid from 380 nm to 780 nm in 5 nm intervals
template<class T> VECTOR<T,3> CIE_XYZ<T>::
Calculate_XYZ(const ARRAY<T,VECTOR<int,1> >& radiometric_radiance_spectrum) const
{
    VECTOR<T,3> XYZ;
    XYZ.y=Radiometric_To_Photometric(radiometric_radiance_spectrum); // Y is the photometric radiance in candela/m^2
    for(int i=1;i<=grid.counts.x;i++){
        XYZ.x+=X_spectrum(i)*radiometric_radiance_spectrum(i);XYZ.z+=Z_spectrum(i)*radiometric_radiance_spectrum(i);}
    T conversion_constant=680,scale=grid.dX.x*conversion_constant;
    XYZ.x*=scale;XYZ.z*=scale;
    return XYZ;
}
//#####################################################################
// Function Calculate_Luminance_Scale_Factors
//#####################################################################
template<class T> void CIE_XYZ<T>:: 
Calculate_Luminance_Scale_Factors()
{
    world_to_display=pow(((T)1.219+pow(display_adaption_luminance,(T).4))/((T)1.219+pow(world_adaption_luminance,(T).4)),(T)2.5);  
    world_to_unit=world_to_display/display_maximum_luminance; 
}
//#####################################################################
// Function World_To_Display_XYZ
//#####################################################################
template<class T> VECTOR<T,3> CIE_XYZ<T>::
World_To_Display_XYZ(const VECTOR<T,3>& XYZ) const
{
    return XYZ*world_to_display;
}
//#####################################################################
// Function World_To_Unit_XYZ
//#####################################################################
template<class T> VECTOR<T,3> CIE_XYZ<T>::
World_To_Unit_XYZ(const VECTOR<T,3>& XYZ) const
{
    return XYZ*world_to_unit;
}
//#####################################################################
// Function XYZ_To_RGB
//#####################################################################
// from http://www.inforamp.net/~poynton/notes/colour_and_gamma/ColorFAQ.html#RTFToC17
// This is for Rec. 709 RGB with the D65 as white point.
template<class T> VECTOR<T,3> CIE_XYZ<T>::   
XYZ_To_RGB(const VECTOR<T,3>& XYZ) const
{
    VECTOR<T,3> RGB;
    RGB.x=(T)3.240479*XYZ.x-(T)1.53715*XYZ.y-(T).498535*XYZ.z; 
    RGB.y=-(T).969256*XYZ.x+(T)1.875991*XYZ.y+(T).041556*XYZ.z; 
    RGB.z=(T).055648*XYZ.x-(T).204043*XYZ.y+(T)1.057311*XYZ.z;
    return RGB;
}
//#####################################################################
// Function RGB_To_XYZ
//#####################################################################   
// from http://www.inforamp.net/~poynton/notes/colour_and_gamma/ColorFAQ.html#RTFToC17
// This is for Rec. 709 RGB with the D65 as white point.
template<class T> VECTOR<T,3> CIE_XYZ<T>::
RGB_To_XYZ(const VECTOR<T,3>& RGB) const
{
    VECTOR<T,3> XYZ;
    XYZ.x=(T).412453*RGB.x+(T).35758*RGB.y+(T).180423*RGB.z;
    XYZ.y=(T).212671*RGB.x+(T).71516*RGB.y+(T).072169*RGB.z;
    XYZ.z=(T).019334*RGB.x+(T).119193*RGB.y+(T).950227*RGB.z;
    return XYZ;
}
//#####################################################################
template class CIE_XYZ<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CIE_XYZ<double>;
#endif
