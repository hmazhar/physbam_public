//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RGB_COLORS 
//#####################################################################
#ifndef __RGB_COLORS__
#define __RGB_COLORS__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class RGB_COLORS
{
public:
    RGB_COLORS() 
    {}

    static VECTOR<T,3> Object_Blue()
    {return VECTOR<double,3>(.3374,.3187,.9);}
    
    static VECTOR<T,3> Ground_Tan()
    {return VECTOR<double,3>(1,.775,.5431);}

    static VECTOR<T,3> Background_Blue()
    {return VECTOR<double,3>(.3137,.3137,.4471);}

    static VECTOR<T,3> Black()
    {return VECTOR<T,3>(0,0,0);}   
    
    static VECTOR<T,3> White()
    {return VECTOR<T,3>(1,1,1);}

    static VECTOR<T,3> Red()
    {return VECTOR<T,3>(1,0,0);}
    
    static VECTOR<T,3> Green()
    {return VECTOR<T,3>(0,1,0);}
    
    static VECTOR<T,3> Blue()
    {return VECTOR<T,3>(0,0,1);}
    
    static VECTOR<T,3> Yellow()
    {return VECTOR<T,3>(1,1,0);}
    
    static VECTOR<T,3> Cyan()
    {return VECTOR<T,3>(0,1,1);}
    
    static VECTOR<T,3> Magenta()
    {return VECTOR<T,3>(1,0,1);}
    
    static VECTOR<T,3> Aquamarine()
    {return VECTOR<double,3>(.439216,.858824,.576471);}
    
    static VECTOR<T,3> Blue_Violet()
    {return VECTOR<double,3>(.62352,.372549,.623529);}
    
    static VECTOR<T,3> Brown()
    {return VECTOR<double,3>(.647059,.164706,.164706);}
    
    static VECTOR<T,3> Cadet_Blue()
    {return VECTOR<double,3>(.372549,.623529,.623529);}
    
    static VECTOR<T,3> Coral()
    {return VECTOR<double,3>(1.0,.498039,0);}
    
    static VECTOR<T,3> Cornflower_Blue()
    {return VECTOR<double,3>(.258824,.258824,.435294);}
    
    static VECTOR<T,3> Dark_Green()
    {return VECTOR<double,3>(.184314,.309804,.184314);}
    
    static VECTOR<T,3> Dark_Olive_Green()
    {return VECTOR<double,3>(.309804,.309804,.184314);}
    
    static VECTOR<T,3> Dark_Orchid()
    {return VECTOR<double,3>(.6,.196078,.8);}
    
    static VECTOR<T,3> Dark_Slate_Blue()
    {return VECTOR<double,3>(.419608,.137255,.556863);}
    
    static VECTOR<T,3> Dark_Slate_Gray()
    {return VECTOR<double,3>(.184314,.309804,.309804);}
    
    static VECTOR<T,3> Dark_Turquoise()
    {return VECTOR<double,3>(.439216,.576471,.858824);}
    
    static VECTOR<T,3> Dim_Gray()
    {return VECTOR<double,3>(.329412,.329412,.329412);}
    
    static VECTOR<T,3> Firebrick()
    {return VECTOR<double,3>(.556863,.137255,.137255);}
    
    static VECTOR<T,3> Forest_Green()
    {return VECTOR<double,3>(.137255,.556863,.137255);}
    
    static VECTOR<T,3> Gold()
    {return VECTOR<double,3>(.8,.498039,.196078);}
    
    static VECTOR<T,3> Goldenrod()
    {return VECTOR<double,3>(.858824,.858824,.439216);}
    
    static VECTOR<T,3> Gray()
    {return VECTOR<double,3>(.752941,.752941,.752941);}
    
    static VECTOR<T,3> Green_Yellow()
    {return VECTOR<double,3>(.576471,.858824,.439216);}
    
    static VECTOR<T,3> Indian_Red()
    {return VECTOR<double,3>(.309804,.184314,.184314);}
    
    static VECTOR<T,3> Khaki()
    {return VECTOR<double,3>(.623529,.623529,.372549);}
    
    static VECTOR<T,3> Light_Blue()
    {return VECTOR<double,3>(.74902,.847059,.847059);}
    
    static VECTOR<T,3> Light_Gray()
    {return VECTOR<double,3>(.658824,.658824,.658824);}
    
    static VECTOR<T,3> Light_Steel_Blue()
    {return VECTOR<double,3>(.560784,.560784,.737255);}
    
    static VECTOR<T,3> Lime_Green()
    {return VECTOR<double,3>(.196078,.8,.196078);}
    
    static VECTOR<T,3> Maroon()
    {return VECTOR<double,3>(.556863,.137255,.419608);}
    
    static VECTOR<T,3> Medium_Aquamarine()
    {return VECTOR<double,3>(.196078,.8,.6);}
    
    static VECTOR<T,3> Medium_Blue()
    {return VECTOR<double,3>(.196078,.196078,.8);}
    
    static VECTOR<T,3> Medium_Forest_Green()
    {return VECTOR<double,3>(.419608,.556863,.137255);}
    
    static VECTOR<T,3> Medium_Goldenrod()
    {return VECTOR<double,3>(.917647,.917647,.678431);}
    
    static VECTOR<T,3> Medium_Orchid()
    {return VECTOR<double,3>(.576471,.439216,.858824);}
    
    static VECTOR<T,3> Medium_Sea_Green()
    {return VECTOR<double,3>(.258824,.435294,.258824);}
    
    static VECTOR<T,3> Medium_Slate_Blue()
    {return VECTOR<double,3>(.498039,0,1);}
    
    static VECTOR<T,3> Medium_Spring_Green()
    {return VECTOR<double,3>(.498039,1.0,0);}
    
    static VECTOR<T,3> Medium_Turquoise()
    {return VECTOR<double,3>(.439216,.858824,.858824);}
    
    static VECTOR<T,3> Medium_Violet_Red()
    {return VECTOR<double,3>(.858824,.439216,.576471);}
    
    static VECTOR<T,3> Midnight_Blue()
    {return VECTOR<double,3>(.184314,.184314,.309804);}
    
    static VECTOR<T,3> Navy()
    {return VECTOR<double,3>(.137255,.137255,.556863);}
    
    static VECTOR<T,3> Navy_Blue()
    {return VECTOR<double,3>(.137255,.137255,.556863);}
    
    static VECTOR<T,3> Orange()
    {return VECTOR<double,3>(.8,.196078,.196078);}
    
    static VECTOR<T,3> Orange_Red()
    {return VECTOR<double,3>(1,0,.498039);}
    
    static VECTOR<T,3> Orchid()
    {return VECTOR<double,3>(.858824,.439216,.858824);}
    
    static VECTOR<T,3> Pale_Green()
    {return VECTOR<double,3>(.560784,.737255,.560784);}
    
    static VECTOR<T,3> Pink()
    {return VECTOR<double,3>(.737255,.560784,.560784);}
    
    static VECTOR<T,3> Plum()
    {return VECTOR<double,3>(.917647,.678431,.917647);}
    
    static VECTOR<T,3> Salmon()
    {return VECTOR<double,3>(.435294,.258824,.258824);}
    
    static VECTOR<T,3> Sea_Green()
    {return VECTOR<double,3>(.137255,.556863,.419608);}
    
    static VECTOR<T,3> Sienna()
    {return VECTOR<double,3>(.556863,.419608,.137255);}
    
    static VECTOR<T,3> Sky_Blue()
    {return VECTOR<double,3>(.196078,.6,.8);}
    
    static VECTOR<T,3> Slate_Blue()
    {return VECTOR<double,3>(0,.498039,1);}
    
    static VECTOR<T,3> Spring_Green()
    {return VECTOR<double,3>(0,1,.498039);}
    
    static VECTOR<T,3> Steel_Blue()
    {return VECTOR<double,3>(.137255,.419608,.556863);}
    
    static VECTOR<T,3> Tan()
    {return VECTOR<double,3>(.858824,.576471,.439216);}
    
    static VECTOR<T,3> Thistle()
    {return VECTOR<double,3>(.847059,.74902,.847059);}
    
    static VECTOR<T,3> Turquoise()
    {return VECTOR<double,3>(.678431,.917647,.917647);}
    
    static VECTOR<T,3> Violet()
    {return VECTOR<double,3>(.309804,.184314,.309804);}
    
    static VECTOR<T,3> Violet_Red()
    {return VECTOR<double,3>(.8,.196078,.6);}
    
    static VECTOR<T,3> Wheat()
    {return VECTOR<double,3>(.847059,.847059,.74902);}
    
    static VECTOR<T,3> Yellow_Green()
    {return VECTOR<double,3>(.6,.8,.196078);}

    static VECTOR<T,3> From_HSV(T h,const T s,const T v)
    {if(s==0) return VECTOR<T,3>(v,v,v);
    h/=360;h-=floor(h);h*=6; // h is now in [0,6)
    int i=(int)floor(h);
    float f=h-i,p=v*(1-s),q=v*(1-s*f),t=v*(1-s*(1-f));
    switch(i){
        case 0:return VECTOR<T,3>(v,t,p);
        case 1:return VECTOR<T,3>(q,v,p);
        case 2:return VECTOR<T,3>(p,v,t);
        case 3:return VECTOR<T,3>(p,q,v);
        case 4:return VECTOR<T,3>(t,p,v);
        case 5:return VECTOR<T,3>(v,p,q);
        default:PHYSBAM_FATAL_ERROR();}}

//#####################################################################
};   
}
#endif
