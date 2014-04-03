//#####################################################################
// Copyright 2002-2003, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BCC_GRID
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __BCC_GRID__
#define __BCC_GRID__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class TV> class RANGE;
template<class T>
class BCC_GRID
{
    typedef VECTOR<T,3> TV;
public:
    T size; // separation distance of primary and secondary grids
    TV base; // origin of bcc grid coordinate system
    T one_over_size,size_over_four,four_over_size;

    BCC_GRID(const T size_input=1)
    {
        Set_Size(size_input);
    }
    
    BCC_GRID(const T size_input,const VECTOR<int,3>& base_input)
        :base(base_input)
    {
        Set_Size(size_input);
    }
    
    void Set_Size(const T size_input)
    {size=size_input;one_over_size=1/size;size_over_four=.25*size;four_over_size=4/size;}
    
    TV Vertex_Position(const VECTOR<int,3>& vertex) const
    {return base+TV(size_over_four*vertex.x,size_over_four*vertex.y,size_over_four*vertex.z);}

    TV Transform_Point(const TV& point) const
    {return four_over_size*(point-base);}
    
    TV Transform_Vector(const TV& vector) const
    {return four_over_size*vector;}
    
    TV Shear_Point(const TV& point) const
    {TV v=one_over_size*(point-base);return TV(v.x+v.y,v.x+v.z,v.y+v.z);}
    
    RANGE<TV> Sheared_Box(const TV& x1,const TV& x2,const TV& x3) const
    {return RANGE<TV>::Bounding_Box(Shear_Point(x1),Shear_Point(x2),Shear_Point(x3));}
    
    VECTOR<int,3> Edge_Between_Vertices(const VECTOR<int,3> v1,const VECTOR<int,3> v2) const
    {return VECTOR<int,3>((v1.x+v2.x)>>1,(v1.y+v2.y)>>1,(v1.z+v2.z)>>1);}
    
    int Level(const VECTOR<int,3>& v) const
    {int o=v.x|v.y|v.z;return o&~(o-1);}
    
    int Maximum_Vertex_Level(const VECTOR<int,3>& v) const
    {int level=Level(v);return (level>>1)-(((v.x^v.y|v.x^v.z)&level)>>2);}
    
    int Edge_Level(const VECTOR<int,3>& edge) const
    {int level=Level(edge);return level-(((edge.x^edge.y|edge.x^edge.z)&level)>>1);}
    
    int Tetrahedron_Level(const VECTOR<int,3>& tet) const
    {return Level(tet);}
    
    bool Valid_Vertex(const VECTOR<int,3>& v,const int level) const
    {int l2=level<<1,mask=l2-1;return !((v.x|v.y|v.z)&mask) && (v.x&l2)==(v.y&l2) && (v.y&l2)==(v.z&l2);}
    
    void Assert_Valid_Edge(const VECTOR<int,3>& e,const int level) const
    {int l2=level<<1,mask=l2-1;assert((e.x&mask)==(e.y&mask) && (e.y&mask)==(e.z&mask) && (e.x&mask || (e.x&l2)!=(e.y&l2) || (e.y&l2)!=(e.z&l2)));}
    
    void Assert_Valid_Tetrahedron(const VECTOR<int,3>& t,const int level) const
    {int l1=level,l2=l1<<1;assert((t.x&l1)+(t.y&l1)+(t.z&l1)==l1 && (t.x&(~t.x<<1)&l2)+(t.y&(~t.y<<1)&l2)+(t.z&(~t.z<<1)&l2)==l2);}
    
    VECTOR<int,3> Tetrahedron_Containing_Point(const TV point,const int level=1) const
    {TV shear=Shear_Point(point);int l=level,mask=~(level-1);
    int ij=(int)shear.x&mask,ik=(int)shear.y&mask,jk=(int)shear.z&mask,i=ij+ik-jk,j=ij+jk-ik,k=ik+jk-ij;
    T xy=shear.x-ij,xz=shear.y-ik,yz=shear.z-jk;int x_gt_y=xz>yz?l:0,x_gt_z=xy>yz?l:0,y_gt_z=xy>xz?l:0;
    return VECTOR<int,3>((i<<1)+x_gt_y+x_gt_z,(j<<1)+l-x_gt_y+y_gt_z,((k+l)<<1)-x_gt_z-y_gt_z);}
    
    void Edges_Intersecting_Sheared_Box(const RANGE<TV>& box,const int level,ARRAY<VECTOR<int,3> >& edge_list) const
    {int ij_min=(int)box.min_corner.x,ij_max=(int)box.max_corner.x,ik_min=(int)box.min_corner.y,ik_max=(int)box.max_corner.y,jk_min=(int)box.min_corner.z,jk_max=(int)box.max_corner.z;
    int l1=level,l2=l1<<1,mask=~(level-1);
    ij_min&=mask;ij_max&=mask;ik_min&=mask;ik_max&=mask;jk_min&=mask;jk_max&=mask;
    for(int ij=ij_min;ij<=ij_max;ij+=l1)for(int ik=ik_min;ik<=ik_max;ik+=l1)for(int jk=jk_min;jk<=jk_max;jk+=l1){
        int i=(ij+ik-jk)<<1,j=(ij+jk-ik)<<1,k=(ik+jk-ij)<<1;
        edge_list.Append(VECTOR<int,3>(i+l2,j,k));edge_list.Append(VECTOR<int,3>(i,j+l2,k));edge_list.Append(VECTOR<int,3>(i,j,k+l2));
        edge_list.Append(VECTOR<int,3>(i+l1,j+l1,k+l1));edge_list.Append(VECTOR<int,3>(i+l1,j+l1,k-l1));
        edge_list.Append(VECTOR<int,3>(i+l1,j-l1,k+l1));edge_list.Append(VECTOR<int,3>(i-l1,j+l1,k+l1));}}
        
    void Tetrahedrons_Intersecting_Sheared_Box(const RANGE<TV>& box,const int level,ARRAY<VECTOR<int,3> >& tetrahedron_list) const
    {int ij_min=(int)box.min_corner.x,ij_max=(int)box.max_corner.x,ik_min=(int)box.min_corner.y,ik_max=(int)box.max_corner.y,jk_min=(int)box.min_corner.z,jk_max=(int)box.max_corner.z;
    int l1=level,l2=l1<<1,mask=~(level-1);
    ij_min&=mask;ij_max&=mask;ik_min&=mask;ik_max&=mask;jk_min&=mask;jk_max&=mask;
    for(int ij=ij_min;ij<=ij_max;ij+=l1)for(int ik=ik_min;ik<=ik_max;ik+=l1)for(int jk=jk_min;jk<=jk_max;jk+=l1){
        int i=(ij+ik-jk)<<1,j=(ij+jk-ik)<<1,k=(ik+jk-ij)<<1;
        tetrahedron_list.Append(VECTOR<int,3>(i+l2,j+l1,k));tetrahedron_list.Append(VECTOR<int,3>(i+l2,j,k+l1));
        tetrahedron_list.Append(VECTOR<int,3>(i+l1,j+l2,k));tetrahedron_list.Append(VECTOR<int,3>(i+l1,j,k+l2));
        tetrahedron_list.Append(VECTOR<int,3>(i,l2+j,l1+k));tetrahedron_list.Append(VECTOR<int,3>(i,l1+j,k+l2));}}
        
    void Vertex_Neighbors(const VECTOR<int,3>& vertex,const int level,ARRAY<VECTOR<int,3> >& vertex_list) const
    {vertex_list.Resize(14);ARRAYS_COMPUTATIONS::Fill(vertex_list,vertex);int l2=level<<1,l4=level<<2;
    vertex_list(1).x+=l4;vertex_list(2).y+=l4;vertex_list(3).z+=l4;vertex_list(4).x-=l4;vertex_list(5).y-=l4;vertex_list(6).z-=l4;
    vertex_list( 7)+=VECTOR<int,3>( l2, l2, l2);vertex_list( 8)+=VECTOR<int,3>( l2, l2,-l2);vertex_list( 9)+=VECTOR<int,3>( l2,-l2, l2);
    vertex_list(10)+=VECTOR<int,3>( l2,-l2,-l2);vertex_list(11)+=VECTOR<int,3>(-l2, l2, l2);vertex_list(12)+=VECTOR<int,3>(-l2, l2,-l2);
    vertex_list(13)+=VECTOR<int,3>(-l2,-l2, l2);vertex_list(14)+=VECTOR<int,3>(-l2,-l2,-l2);}
    
    void Incident_Tetrahedrons(const VECTOR<int,3>& vertex,const int level,ARRAY<VECTOR<int,3> >& tetrahedron_list) const
    {ARRAY<VECTOR<int,3> >& tl=tetrahedron_list;tl.Resize(24);ARRAYS_COMPUTATIONS::Fill(tl,vertex);int t=1;
    int l1=level,l2=l1<<1;
    tl( 1).x+=l2;tl( 1).y+=l1;tl( 2).x+=l2;tl( 2).y-=l1;tl( 3).x+=l2;tl( 3).z+=l1;tl( 4).x+=l2;tl( 4).z-=l1;
    tl( 5).x-=l2;tl( 5).y+=l1;tl( 6).x-=l2;tl( 6).y-=l1;tl( 7).x-=l2;tl( 7).z+=l1;tl( 8).x-=l2;tl( 8).z-=l1;
    tl( 9).y+=l2;tl( 9).x+=l1;tl(10).y+=l2;tl(10).x-=l1;tl(11).y+=l2;tl(11).z+=l1;tl(12).y+=l2;tl(12).z-=l1;
    tl(13).y-=l2;tl(13).x+=l1;tl(14).y-=l2;tl(14).x-=l1;tl(15).y-=l2;tl(15).z+=l1;tl(16).y-=l2;tl(16).z-=l1;
    tl(17).z+=l2;tl(17).x+=l1;tl(18).z+=l2;tl(18).x-=l1;tl(19).z+=l2;tl(19).y+=l1;tl(20).z+=l2;tl(20).y-=l1;
    tl(21).z-=l2;tl(21).x+=l1;tl(22).z-=l2;tl(22).x-=l1;tl(23).z-=l2;tl(23).y+=l1;tl(24).z-=l2;tl(24).y-=l1;}
        
    void Edge_Vertices(const VECTOR<int,3>& edge,VECTOR<int,3>& v1,VECTOR<int,3>& v2) const
    {int l1=Edge_Level(edge),l2=l1<<1,xp=edge.x&l2,yp=edge.y&l2,zp=edge.z&l2;VECTOR<int,3> d;
    d=edge.x&l1?VECTOR<int,3>(xp-l1,yp-l1,zp-l1):VECTOR<int,3>(yp^zp^l2,xp^zp^l2,xp^yp^l2);v1=edge-d;v2=edge+d;}
    
    void Edge_Tetrahedrons(const VECTOR<int,3>& edge,ARRAY<VECTOR<int,3> >& tetrahedron_list) const
    {int l1=Edge_Level(edge),l2=l1<<1;Assert_Valid_Edge(edge,l1);
    int xp=edge.x&l2,yp=edge.y&l2,zp=edge.z&l2,t=1;
    if(edge.x&l1){
        tetrahedron_list.Resize(6);ARRAYS_COMPUTATIONS::Fill(tetrahedron_list,edge);
        int dx=l1-xp,dy=l1-yp,dz=l1-zp;
        tetrahedron_list(1).x+=dx;tetrahedron_list(1).y-=dy;tetrahedron_list(2).x+=dx;tetrahedron_list(2).z-=dz;
        tetrahedron_list(3).y+=dy;tetrahedron_list(3).x-=dx;tetrahedron_list(4).y+=dy;tetrahedron_list(4).z-=dz;
        tetrahedron_list(5).z+=dz;tetrahedron_list(5).x-=dx;tetrahedron_list(6).z+=dz;tetrahedron_list(6).y-=dy;}
    else{
        tetrahedron_list.Resize(4);ARRAYS_COMPUTATIONS::Fill(tetrahedron_list,edge);
        if(yp^zp){tetrahedron_list(t++).x+=l1;tetrahedron_list(t++).x-=l1;}
        if(xp^zp){tetrahedron_list(t++).y+=l1;tetrahedron_list(t++).y-=l1;}
        if(xp^yp){tetrahedron_list(t++).z+=l1;tetrahedron_list(t++).z-=l1;}
        assert(t==5);}
    for(int i=1;i<=tetrahedron_list.m;i++)Assert_Valid_Tetrahedron(tetrahedron_list(i),l1);}
    
    void Edge_Children(const VECTOR<int,3>& edge,VECTOR<int,3>& e1,VECTOR<int,3>& e2) const
    {int l2=Edge_Level(edge),l1=l2>>1,xp=(edge.x>>1)&l2,yp=(edge.y>>1)&l2,zp=(edge.z>>1)&l2;VECTOR<int,3> d;
    d=edge.x&l2?VECTOR<int,3>(xp-l1,yp-l1,zp-l1):VECTOR<int,3>(yp^zp^l2,xp^zp^l2,xp^yp^l2);e1=edge-d;e2=edge+d;}
    
    int Edge_Parent(const VECTOR<int,3>& edge,VECTOR<int,3>& parent) const
    {int l1=Edge_Level(edge),l2=l1<<1,l4=l1<<2,x2=edge.x&l2,y2=edge.y&l2,z2=edge.z&l2,x4=edge.x&l4,y4=edge.y&l4,z4=edge.z&l4;
    if(edge.x&l1){
        int xs=x4^(x2<<1);
        if(xs^y4^(y2<<1) || xs^z4^(z2<<1))return 0;
        parent=edge+VECTOR<int,3>(l1-x2,l1-y2,l1-z2);return 1-((xs!=0)<<1);}
    else{parent=edge;
        if(x2+y2+z2!=l2)return 0;
        if(x2){if(y4!=z4)return 0;parent.x+=l2-(x4^y4);return 1;}
        else if(y2){if(z4!=x4)return 0;parent.y+=l2-(y4^z4);return 1;}
        else{assert(z2);if(x4!=y4)return 0;parent.z+=l2-(z4^x4);return 1;}}}
        
    void Hierarchy_Check(VECTOR<int,3>& edge) const
    {int level=Edge_Level(edge);Assert_Valid_Edge(edge,level);
    VECTOR<int,3> parent;int s=Edge_Parent(edge,parent);
    VECTOR<int,3> e1,e2;
    if(s && s < 4){
        VECTOR<int,3> parent_v1,parent_v2,edge_v1,edge_v2;
        Assert_Valid_Edge(parent,level<<1);assert(-1<=s && s<=1);
        Edge_Children(parent,e1,e2);
        assert(edge==e1 || edge==e2);
        Edge_Vertices(parent,parent_v1,parent_v2);Edge_Vertices(edge,edge_v1,edge_v2);
        assert(edge_v1==parent||edge_v2==parent);
        if(s>0) assert(parent_v1==edge_v1||parent_v2==edge_v2);
        else assert(parent_v1==edge_v2||parent_v2==edge_v1);}
    if(level>1){
        Edge_Children(edge,e1,e2);
        assert(Edge_Parent(e1,parent) && Edge_Parent(e2,parent));}}

    void Tetrahedron_Vertices(const VECTOR<int,3>& t,VECTOR<int,3>& v1,VECTOR<int,3>& v2,VECTOR<int,3>& v3,VECTOR<int,3>& v4) const
    {int l1=Tetrahedron_Level(t),l2=l1<<1;
    if(t.x&l1){
        int s=l1-((t.x^t.y)&l2),s2=s<<1;
        v1=t+VECTOR<int,3>(s,l2,0);v2=t+VECTOR<int,3>(s,-l2,0);
        v3=t+VECTOR<int,3>(-s,0,s2);v4=t+VECTOR<int,3>(-s,0,-s2);}
    else if(t.y&l1){
        int s=l1-((t.y^t.z)&l2),s2=s<<1;
        v1=t+VECTOR<int,3>(0,s,l2);v2=t+VECTOR<int,3>(0,s,-l2);
        v3=t+VECTOR<int,3>(s2,-s,0);v4=t+VECTOR<int,3>(-s2,-s,0);}
    else{assert(t.z&l1);
        int s=l1-((t.z^t.x)&l2),s2=s<<1;
        v1=t+VECTOR<int,3>(l2,0,s);v2=t+VECTOR<int,3>(-l2,0,s);
        v3=t+VECTOR<int,3>(0,s2,-s);v4=t+VECTOR<int,3>(0,-s2,-s);}}

    void Tetrahedron_Children(const VECTOR<int,3>& t,ARRAY<VECTOR<int,3> >& child_list) const
    {ARRAY<VECTOR<int,3> >&cl=child_list;cl.Resize(8);ARRAYS_COMPUTATIONS::Fill(cl,t);
    int l2=Tetrahedron_Level(t),l4=l2<<1,l1=l2>>1;
    if(t.x&l2){
        int s=(l2-((t.x^t.y)&l4))>>1;
        cl(1).y-=l1;cl(2).y+=l1;cl(3).z-=l1;cl(4).z+=l1;
        cl(5).x+=s;cl(5).y-=l2;cl(6).x+=s;cl(6).y+=l2;cl(7).x-=s;cl(7).z-=l2;cl(8).x-=s;cl(8).z+=l2;}
    else if(t.y&l2){
        int s=(l2-((t.y^t.z)&l4))>>1;
        cl(1).z-=l1;cl(2).z+=l1;cl(3).x-=l1;cl(4).x+=l1;
        cl(5).y+=s;cl(5).z-=l2;cl(6).y+=s;cl(6).z+=l2;cl(7).y-=s;cl(7).x-=l2;cl(8).y-=s;cl(8).x+=l2;}
    else{assert(t.z&l2);
        int s=(l2-((t.z^t.x)&l4))>>1;
        cl(1).x-=l1;cl(2).x+=l1;cl(3).y-=l1;cl(4).y+=l1;
        cl(5).z+=s;cl(5).x-=l2;cl(6).z+=s;cl(6).x+=l2;cl(7).z-=s;cl(7).y-=l2;cl(8).z-=s;cl(8).y+=l2;}}
        
    void Tetrahedron_Parent(const VECTOR<int,3>& t,VECTOR<int,3>& parent) const
    {int l1=Tetrahedron_Level(t),l2=l1<<1,l4=l2<<1;parent=t;
    if(t.x&l1){
        int xp=t.x&l2,xp2=xp<<1,s=l1-xp,ss=l2-((t.y^t.z)&l4);
        if(t.y&l2){
            if((xp2^t.x^t.z)&l4){parent.x-=s;}
            else{parent.x+=s;parent.y+=ss;}}
        else{
            if((xp2^t.x^t.y)&l4){parent.x-=s;}
            else{parent.x+=s;parent.z+=ss;}}}
    else if(t.y&l1){
        int yp=t.y&l2,yp2=yp<<1,s=l1-yp,ss=l2-((t.z^t.x)&l4);
        if(t.z&l2){
            if((yp2^t.y^t.x)&l4){parent.y-=s;}
            else{parent.y+=s;parent.z+=ss;}}
        else{
            if((yp2^t.y^t.z)&l4){parent.y-=s;}
            else{parent.y+=s;parent.x+=ss;}}}
    else{assert(t.z&l1);
        int zp=t.z&l2,zp2=zp<<1,s=l1-zp,ss=l2-((t.x^t.y)&l4);
        if(t.x&l2){
            if((zp2^t.z^t.y)&l4){parent.z-=s;}
            else{parent.z+=s;parent.x+=ss;}}
        else{
            if((zp2^t.z^t.x)&l4){parent.z-=s;}
            else{parent.z+=s;parent.y+=ss;}}}}
            
//#####################################################################
};
}
#endif
#endif
