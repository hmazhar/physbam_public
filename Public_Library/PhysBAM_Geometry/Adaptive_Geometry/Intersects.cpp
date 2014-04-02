//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_FLOAT.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Intersects.h>
#include <PhysBAM_Geometry/Fracture/CUTTING_SIMPLEX.h>
#include <PhysBAM_Geometry/Fracture/CUTTING_SIMPLEX_COORDINATE.h>
namespace PhysBAM{
namespace GEOMETRIC_PREDICATES_DETAIL{
    //#####################################################################
    // Has_Separating_Hyperplane
    //#####################################################################
    template<class T_EXACT,class T,int m,int n> bool
    Has_Separating_Hyperplane(const VECTOR<VECTOR<T,1>,m>& simplex1,const VECTOR<VECTOR<T,1>,n>& simplex2,bool* is_degenerate_p)
    {
        PHYSBAM_STATIC_ASSERT((m<=2&&n<=2));
        typedef VECTOR<T,1> TV;
        if(is_degenerate_p) *is_degenerate_p=false;
        bool is_possibly_degenerate=false;
        // attempt to find a separator point composed of 1 vertex from the union
        //  of vertices of both simplices
        VECTOR<TV,m+n> all_vertices=simplex1.Append_Elements(simplex2);
        for(int v1=1;v1<m+n;++v1){
            const TV& x1=all_vertices[v1];
            // x1 defines the trial separator point
            // sign[i] indicates which side of the trial point simplexi's vertices
            //  fall on
            bool could_be_a_separator=true;
            VECTOR<int,2> sign;
            int next_sign=0;
            for(int u=1;u<=m+n&&could_be_a_separator;++u){
                if(u==v1) continue;
                const TV& y=all_vertices[u];
                next_sign=Adaptive_Signed_Volume<T_EXACT,T>(x1,y).Sign();
                int i=(u<=m?1:2);
                if(sign[i]==0) sign[i]=next_sign;
                could_be_a_separator=(next_sign!=0&&sign[i]==next_sign);}
            if(could_be_a_separator){
                assert(sign[1]!=0||sign[2]!=0);
                if(sign[1]!=sign[2]) return true;}
            else if(next_sign==0) is_possibly_degenerate=true;}
        if(is_possibly_degenerate){
            if(is_degenerate_p==0) throw GEOMETRIC_DEGENERACY();
            *is_degenerate_p=true;}
        return false;
    }    
    template<class T_EXACT,class T,int m,int n> bool
    Has_Separating_Hyperplane(const VECTOR<VECTOR<T,2>,m>& simplex1,const VECTOR<VECTOR<T,2>,n>& simplex2,bool* is_degenerate_p)
    {
        PHYSBAM_STATIC_ASSERT((m<=3&&n<=3));
        typedef VECTOR<T,2> TV;
        if(is_degenerate_p) *is_degenerate_p=false;
        bool is_possibly_degenerate=false;
        // attempt to find a separator line composed of 2 vertices from the union
        //  of vertices of both simplices
        VECTOR<TV,m+n> all_vertices=simplex1.Append_Elements(simplex2);
        for(int v1=1;v1<m+n;++v1){
            const TV& x1=all_vertices[v1];
            for(int v2=v1+1;v2<=m+n;++v2){
                const TV& x2=all_vertices[v2];
                // x1,x2 defines the trial separator line
                // sign[i] indicates which side of the trial line simplexi's
                //  vertices fall on
                bool could_be_a_separator=true;
                VECTOR<int,2> sign;
                int next_sign=0;
                for(int u=1;u<=m+n&&could_be_a_separator;++u){
                    if(u==v1||u==v2) continue;
                    const TV& y=all_vertices[u];
                    next_sign=Adaptive_Signed_Volume<T_EXACT>(x1,x2,y).Sign();
                    int i=(u<=m?1:2);
                    if(sign[i]==0) sign[i]=next_sign;
                    could_be_a_separator=(next_sign!=0&&sign[i]==next_sign);}
                if(could_be_a_separator){
                    assert(sign[1]!=0||sign[2]!=0);
                    if(sign[1]!=sign[2]) return true;}
                else if(next_sign==0) is_possibly_degenerate=true;}}
        if(is_possibly_degenerate){
            if(is_degenerate_p==0) throw GEOMETRIC_DEGENERACY();
            *is_degenerate_p=true;}
        return false;
    }
    
    template<class T_EXACT,class T,int m,int n> bool
    Has_Separating_Hyperplane(const VECTOR<VECTOR<T,3>,m>& simplex1,const VECTOR<VECTOR<T,3>,n>& simplex2,bool* is_degenerate_p)
    {
        PHYSBAM_STATIC_ASSERT((m<=4&&n<=4));
        typedef VECTOR<T,3> TV;
        if(is_degenerate_p) *is_degenerate_p=false;
        bool is_possibly_degenerate=false;
        // attempt to find a separator plane composed of 3 vertices from the union
        //  of vertices of both simplices
        VECTOR<TV,m+n> all_vertices=simplex1.Append_Elements(simplex2);
        for(int v1=1;v1<m+n;++v1){
            const TV& x1=all_vertices[v1];
            for(int v2=v1+1;v2<m+n;++v2){
                const TV& x2=all_vertices[v2];
                for(int v3=v2+1;v3<=m+n;++v3){
                    const TV& x3=all_vertices[v3];
                    // x1,x2,x3 defines the trial separator plane
                    // sign[i] indicates which side of the trial plane simplexi's
                    //  vertices fall on
                    bool could_be_a_separator=true;
                    VECTOR<int,2> sign(0,0);
                    int next_sign=0;
                    for(int u=1;u<=m+n&&could_be_a_separator;++u){
                        if(u==v1||u==v2||u==v3) continue;
                        const TV& y=all_vertices[u];
                        next_sign=Adaptive_Signed_Volume<T_EXACT>(x1,x2,x3,y).Sign();
                        int i=(u<=m?1:2);
                        if(sign[i]==0) sign[i]=next_sign;
                        could_be_a_separator=(next_sign!=0&&sign[i]==next_sign);}
                    if(could_be_a_separator){
                        assert(sign[1]!=0||sign[2]!=0);
                        if(sign[1]!=sign[2]) return true;}
                    else if(next_sign==0) is_possibly_degenerate=true;}}}
        if(is_possibly_degenerate){
            if(is_degenerate_p==0) throw GEOMETRIC_DEGENERACY();
            *is_degenerate_p=true;}
        return false;
    }

    //#####################################################################
    // Is_Cohyperplanar
    //#####################################################################
    template<class T_EXACT,class T,int d,int m,int n> bool
    Is_Cohyperplanar(const VECTOR<VECTOR<T,d>,m>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2)
    {
        PHYSBAM_STATIC_ASSERT((m<=d&&n<=d));
        typedef VECTOR<T,d> TV;
        VECTOR<TV,m+n> all_vertices=simplex1.Append_Elements(simplex2);
        VECTOR<TV,d> matrix;
        for(int i=1;i<=d;++i) matrix[i]=all_vertices[i];
        for(int i=d+1;i<=m+n;++i) if(Adaptive_Signed_Volume<T_EXACT>(all_vertices[i],matrix).Sign()!=0) return false;
        return true;
    }
    
    template<class T_EXACT,class T,int d,int m,int n> bool
    Intersects_Or_Is_Degenerate_After_Projection(const VECTOR<VECTOR<T,d>,m>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2)
    {
        PHYSBAM_STATIC_ASSERT((m<=d&&n<=d));
        assert(Is_Cohyperplanar<T_EXACT>(simplex1,simplex2));
        VECTOR<VECTOR<T,d-1>,m> simplex1_proj;
        VECTOR<VECTOR<T,d-1>,n> simplex2_proj;
        bool intersects=false,is_degenerate=true;
        for(int k=1;k<=d&&is_degenerate;++k){
            for(int i=1;i<=m;++i) simplex1_proj[i]=simplex1[i].Remove_Index(k);
            for(int i=1;i<=n;++i) simplex2_proj[i]=simplex2[i].Remove_Index(k);
            intersects=Intersects<T_EXACT>(simplex1_proj,simplex2_proj,&is_degenerate);}
        return (intersects||is_degenerate);
    }

    //#####################################################################
    // Is_Any_Inside
    //#####################################################################
    template<class T_EXACT,class T,int d,int n> bool
    Is_Any_Inside_Helper(const VECTOR<VECTOR<T,d>,d+1>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2)
    {bool is_degenerate;
    for(int i=1;i<=n;++i) if(Intersects<T_EXACT>(simplex1,simplex2[i],&is_degenerate)&&!is_degenerate) return true;
    return false;}

template<class T_EXACT,class T,int d> bool
Intersects(const VECTOR<VECTOR<T,d>,d+1>& simplex,const VECTOR<T,d>& point,bool* is_degenerate_p)
{
    VECTOR<VECTOR<T,d>,d+1> test_simplex=simplex;
    int sign=0;
    bool is_possibly_degenerate=false;
    if(is_degenerate_p) *is_degenerate_p=false;
    // compute the signed volumes obtained by replacing each simplex vertex in turn with point if 2 such signed volumes have opposite sign, point cannot lie in simplex
    for(int i=1;i<=d+1;++i){
        test_simplex[i]=point;
        int next_sign=Adaptive_Signed_Volume<T_EXACT>(test_simplex).Sign();
        if(next_sign==0) is_possibly_degenerate=true;
        else if(sign==0) sign=next_sign;
        else if(next_sign!=sign) return false;
        test_simplex[i]=simplex[i];}
    if(is_possibly_degenerate){if(is_degenerate_p==0) throw GEOMETRIC_DEGENERACY();*is_degenerate_p=true;}
    return true;
}
template<class T_EXACT,class T,int d,int m,int n> typename ENABLE_IF<(n<=m&&m<=d&&m+n<=d+2),bool>::TYPE
Intersects(const VECTOR<VECTOR<T,d>,m>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2,bool* is_degenerate_p)
{
    if(is_degenerate_p) *is_degenerate_p=false;
    bool is_degenerate;
    // check for "normal" intersection first
    bool intersects=!Has_Separating_Hyperplane<T_EXACT>(simplex1,simplex2,&is_degenerate);
    if(!is_degenerate) return intersects;
    // check if coplanar, and if so, perform intersection check on one dimension lower
    if(Is_Cohyperplanar<T_EXACT>(simplex1,simplex2)&&!Intersects_Or_Is_Degenerate_After_Projection<T_EXACT>(simplex1,simplex2)) return false;
    // do not pursue this further, for now
    if(is_degenerate_p==0) throw GEOMETRIC_DEGENERACY();
    *is_degenerate_p=true;
    return true;
}

template<class T_EXACT,class T,int d,int m,int n> typename ENABLE_IF<(n<=m&&m<=d&&m+n>=d+3),bool>::TYPE
Intersects(const VECTOR<VECTOR<T,d>,m>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2,bool* is_degenerate_p)
{
    if(is_degenerate_p!=0) *is_degenerate_p=false;
    bool is_degenerate;
    // check for "normal" intersection first
    bool intersects=!Has_Separating_Hyperplane<T_EXACT>(simplex1,simplex2,&is_degenerate);
    if(!is_degenerate) return intersects;
    // check if coplanar, and if so, perform intersection check on one dimension lower
    if(Is_Cohyperplanar<T_EXACT>(simplex1,simplex2)){if(!Intersects_Or_Is_Degenerate_After_Projection<T_EXACT>(simplex1,simplex2)) return false;}
    else for(int i=1;i<=m;++i) if(Intersects<T_EXACT>(simplex1.Remove_Index(i),simplex2,&is_degenerate)&&!is_degenerate) return true; // check if a subsimplex of simplex1 intersects simplex2
    // do not pursue this further, for now
    if(is_degenerate_p==0) throw GEOMETRIC_DEGENERACY();
    *is_degenerate_p=true;
    return true;
}

template<class T_EXACT,class T,int d,int n> bool
Intersects(const VECTOR<VECTOR<T,d>,d+1>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2,bool* is_degenerate_p)
{
    if(is_degenerate_p) *is_degenerate_p=false;
    bool is_degenerate;
    // check for "normal" intersection first
    bool intersects=!Has_Separating_Hyperplane<T_EXACT>(simplex1,simplex2,&is_degenerate);
    if(!is_degenerate) return intersects;
    // check if a point from either simplex lies strictly inside the other
    if(Is_Any_Inside<T_EXACT>(simplex1,simplex2)) return true;
    // check if a subsimplex of simplex1 intersects simplex2
    for(int i=1;i<=d+1;++i) if(Intersects<T_EXACT>(simplex1.Remove_Index(i),simplex2,&is_degenerate)&&!is_degenerate) return true;
    // do not pursue this further, for now
    if(is_degenerate_p==0) throw GEOMETRIC_DEGENERACY();
    *is_degenerate_p=true;
    return true;
}

template<class T_EXACT,class T>
bool Intersects(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,bool* is_degenerate_p)
{
    typedef typename INTERSECTION_COORDINATES_ADAPTIVE_RESULT<T_EXACT,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3>,VECTOR<VECTOR<T,3>,3> >::TYPE T_ADAPTIVE;
    if(is_degenerate_p!=0) *is_degenerate_p=false;
    VECTOR<T_ADAPTIVE,3> triangle1_coordinates,triangle2_coordinates,triangle3_coordinates;
    return Intersects_And_Intersection_Coordinates<T_EXACT>(triangle1,triangle2,triangle3,triangle1_coordinates,triangle2_coordinates,triangle3_coordinates,is_degenerate_p);
}

template ENABLE_IF<(((2) <= (2)) && ((2) <= (2))) && (((2) + (2)) <= ((2) + (2))),bool>::TYPE Intersects<EXACT_FLOAT<float>,float,2,2,2>(VECTOR<VECTOR<float,2>,2> const&,VECTOR<VECTOR<float,2>,2> const&,bool*);
template ENABLE_IF<(((2) <= (3)) && ((3) <= (3))) && (((3) + (2)) <= ((3) + (2))),bool>::TYPE Intersects<void,CUTTING_SIMPLEX_COORDINATE<float,3>,3,3,2>(VECTOR<VECTOR<CUTTING_SIMPLEX_COORDINATE<float,3>,3>,3> const&,VECTOR<VECTOR<CUTTING_SIMPLEX_COORDINATE<float,3>,3>,2> const&,bool*);
template ENABLE_IF<(((3) <= (3)) && ((3) <= (3))) && (((3) + (3)) >= ((3) + (3))),bool>::TYPE Intersects<void,CUTTING_SIMPLEX_COORDINATE<float,3>,3,3,3>(VECTOR<VECTOR<CUTTING_SIMPLEX_COORDINATE<float,3>,3>,3> const&,VECTOR<VECTOR<CUTTING_SIMPLEX_COORDINATE<float,3>,3>,3> const&,bool*);
template bool Intersects<EXACT_FLOAT<float>,float,2>(VECTOR<VECTOR<float,2>,(2) + (1)> const&,VECTOR<float,2> const&,bool*);
template bool Intersects<EXACT_FLOAT<float>,float,3,2>(VECTOR<VECTOR<float,3>,(3) + (1)> const&,VECTOR<VECTOR<float,3>,2> const&,bool*);
template bool Intersects<EXACT_FLOAT<float>,float,3,3>(VECTOR<VECTOR<float,3>,(3) + (1)> const&,VECTOR<VECTOR<float,3>,3> const&,bool*);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template ENABLE_IF<(((2) <= (2)) && ((2) <= (2))) && (((2) + (2)) <= ((2) + (2))),bool>::TYPE Intersects<EXACT_FLOAT<double>,double,2,2,2>(VECTOR<VECTOR<double,2>,2> const&,VECTOR<VECTOR<double,2>,2> const&,bool*);
template ENABLE_IF<(((2) <= (3)) && ((3) <= (3))) && (((3) + (2)) <= ((3) + (2))),bool>::TYPE Intersects<void,CUTTING_SIMPLEX_COORDINATE<double,3>,3,3,2>(VECTOR<VECTOR<CUTTING_SIMPLEX_COORDINATE<double,3>,3>,3> const&,VECTOR<VECTOR<CUTTING_SIMPLEX_COORDINATE<double,3>,3>,2> const&,bool*);
template ENABLE_IF<(((3) <= (3)) && ((3) <= (3))) && (((3) + (3)) >= ((3) + (3))),bool>::TYPE Intersects<void,CUTTING_SIMPLEX_COORDINATE<double,3>,3,3,3>(VECTOR<VECTOR<CUTTING_SIMPLEX_COORDINATE<double,3>,3>,3> const&,VECTOR<VECTOR<CUTTING_SIMPLEX_COORDINATE<double,3>,3>,3> const&,bool*);
template bool Intersects<EXACT_FLOAT<double>,double,2>(VECTOR<VECTOR<double,2>,(2) + (1)> const&,VECTOR<double,2> const&,bool*);
template bool Intersects<EXACT_FLOAT<double>,double,3,2>(VECTOR<VECTOR<double,3>,(3) + (1)> const&,VECTOR<VECTOR<double,3>,2> const&,bool*);
template bool Intersects<EXACT_FLOAT<double>,double,3,3>(VECTOR<VECTOR<double,3>,(3) + (1)> const&,VECTOR<VECTOR<double,3>,3> const&,bool*);
#endif
}
}
