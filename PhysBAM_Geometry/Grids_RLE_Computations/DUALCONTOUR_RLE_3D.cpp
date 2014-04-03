//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DUALCONTOUR_RLE_3D
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Grids_RLE_Computations/DUALCONTOUR_RLE_3D.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Function Dualcontour
//#####################################################################
template<class T> void DUALCONTOUR_RLE_3D<T>::
Dualcontour()
{
    levelset.Compute_Normals();
    // generate quadrilateral indices (negative indices represent reverse orientation)
    int total_quads=0;
    ARRAY<int> face_quads(levelset.grid.number_of_faces);
    RLE_GRID_3D<T>::template Face_Loop<Generate_Quadrilateral_Indices>(*this,face_quads,total_quads);
    quads.Exact_Resize(total_quads);
    // generate vertices and fill in quad vertices
    positions_and_normals.Preallocate(total_quads);
    for(RLE_GRID_ITERATOR_BLOCK_3D<T> block(levelset.grid,number_of_ghost_cells);block;block++){
        int fx[4];fx[0]=block.Face_X(1);fx[1]=fx[0]+1;fx[3]=block.Face_X(7);fx[2]=fx[3]+1; // counterclockwise around x
        int fz[4];fz[0]=block.Face_Z(4);fz[1]=block.Face_Z(5);fz[3]=fz[0]+1;fz[2]=fz[1]+1; // counterclockwise around z
        int fy[4];fy[0]=block.Face_Y(0)+1;fy[1]=block.Face_Y(6)+1;fy[2]=block.Face_Y(7)+1;fy[3]=block.Face_Y(1)+1; // counterclockwise around y
        int q[3][4];bool need_vertex=false;
        for(int i=0;i<4;i++)if((q[0][i]=face_quads(fx[i]))) need_vertex=true;
        for(int i=0;i<4;i++)if((q[1][i]=face_quads(fy[i]))) need_vertex=true;
        for(int i=0;i<4;i++)if((q[2][i]=face_quads(fz[i]))) need_vertex=true;
        if(!need_vertex) continue;
        int vertex=positions_and_normals.m+1;
        // fill in vertex in quad accounting for orientation of quad
        for(int d=0;d<3;d++)for(int i=0;i<4;i++)if(q[d][i]){
            if(q[d][i]>0) quads(q[d][i])(((i+2)&3)+1)=vertex;
            else quads(-q[d][i])(((5-i)&3)+1)=vertex;}
        // generate vertex
        TV center=block.Center(),position=center;
        TV normal=levelset.Normal(block,position);
        T phi=levelset.Phi(block,position);
        int iterations=0;
        while(abs(phi)>1e-6 && (iterations++)<10){
            position-=phi*normal;
            RLE_GRID_ITERATOR_BLOCK_3D<T> new_block(levelset.grid,position);
            while(!new_block){position=(T).5*(position+center);new_block.Initialize(position);}
            phi=levelset.Phi(new_block,position);
            normal=levelset.Normal(new_block,position);}
        // project vertex back into block, TODO: replace this with something much, much smarter
        BOX<TV> box(block.Bounding_Box());
        box.Scale_About_Center(2); // removes some of the artifacts caused by the following line
        if(box.Lazy_Outside(position)) position=box.Surface(position);
        positions_and_normals.Append(VECTOR<TV,2>(position,normal));}
    // remove incomplete quads
    for(int q=quads.m;q>=1;q--){
        int i,j,k,l;quads(q).Get(i,j,k,l);
        if(!(i&&j&&k&&l)) quads.Remove_Index_Lazy(q);}
}
//#####################################################################
// Function Generate_Quadrilateral_Indices
//#####################################################################
template<class T> template<class T_FACE> void DUALCONTOUR_RLE_3D<T>::
Generate_Quadrilateral_Indices::Apply(const DUALCONTOUR_RLE_3D<T>& dualcontour,ARRAY<int>& face_quads,int& total_quads)
{
    const ARRAY<T>& phi=dualcontour.levelset.phi;
    for(T_FACE face(dualcontour.levelset.grid,dualcontour.number_of_ghost_cells,false);face;face++)if(face.Both_Cells_Short()){int f=face.Face(),c1=face.Cell(0),c2=face.Cell(1);
        if(LEVELSET_UTILITIES<T>::Interface(phi(c1),phi(c2))){total_quads++;face_quads(f)=phi(c2)>0?total_quads:-total_quads;}}
}
//#####################################################################
// Function Get_Triangulated_Surface
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* DUALCONTOUR_RLE_3D<T>::
Get_Triangulated_Surface()
{
    TRIANGULATED_SURFACE<T>* triangulated_surface=TRIANGULATED_SURFACE<T>::Create();
    TRIANGLE_MESH& mesh=triangulated_surface->mesh;
    mesh.number_nodes=positions_and_normals.m;
    mesh.elements.Exact_Resize(2*quads.m);
    int current_triangle=1;
    for(int q=1;q<=quads.m;q++){
        int i,j,k,l;quads(q).Get(i,j,k,l);
        mesh.elements(current_triangle++).Set(i,j,l);mesh.elements(current_triangle++).Set(j,k,l);}
    GEOMETRY_PARTICLES<TV>& particles=triangulated_surface->particles;
    particles.array_collection->Add_Elements(positions_and_normals.m);
    ARRAY<TV>* vertex_normals=new ARRAY<TV>(positions_and_normals.m);
    for(int p=1;p<=positions_and_normals.m;p++) positions_and_normals(p).Get(particles.X(p),(*vertex_normals)(p));
    triangulated_surface->Use_Vertex_Normals();triangulated_surface->vertex_normals=vertex_normals;
    return triangulated_surface;
}
//#####################################################################
template class DUALCONTOUR_RLE_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DUALCONTOUR_RLE_3D<double>;
#endif
#endif
