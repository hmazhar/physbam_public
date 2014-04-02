//#####################################################################
// Copyright 2004-2008, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINTS_POOL
//#####################################################################
#ifndef __POINTS_POOL__
#define __POINTS_POOL__

#include <PhysBAM_Tools/Clone/CLONE_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Log/LOG.h>
namespace PhysBAM{

template<class T_POINT_CLOUD>
class POINTS_POOL:public NONCOPYABLE
{
private:
    typedef typename T_POINT_CLOUD::VECTOR_T TV;

    const T_POINT_CLOUD& template_particles;
    ARRAY<CLONE_ARRAY<T_POINT_CLOUD>*> allocated_batches; // pointer to the first particle in the batch (of size allocation_batch_size)
    STACK<T_POINT_CLOUD*> free_pool;
    int allocation_batch_size;
#ifdef USE_PTHREADS
    pthread_mutex_t stack_lock,batch_lock;
#endif
public:
    int number_particles_per_cell;

    POINTS_POOL(const T_POINT_CLOUD& template_particles,const int allocation_batch_size_input=100000)
        :template_particles(template_particles),allocation_batch_size(allocation_batch_size_input),number_particles_per_cell(0)
    {
#ifdef USE_PTHREADS
        pthread_mutex_init(&stack_lock,0);pthread_mutex_init(&batch_lock,0);
#endif
    }

    ~POINTS_POOL()
    {allocated_batches.Delete_Pointers_And_Clean_Memory();}

    void Set_Number_Particles_Per_Cell(const int number_particles_per_cell_input)
    {number_particles_per_cell=number_particles_per_cell_input;}

    T_POINT_CLOUD* Allocate_Particle()
    {if(free_pool.Empty()){
#ifdef USE_PTHREADS
        pthread_mutex_lock(&batch_lock);
        if(free_pool.Empty()) Allocate_New_Batch();
        pthread_mutex_unlock(&batch_lock);
#else
        Allocate_New_Batch();
#endif
    }
#ifdef USE_PTHREADS
    pthread_mutex_lock(&stack_lock);
    if(free_pool.Empty()){pthread_mutex_unlock(&stack_lock);return Allocate_Particle();} //Some other thread got what was there so try again
#endif
    T_POINT_CLOUD* particles=free_pool.Pop();
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&stack_lock);
#endif
    return particles;}

    void Free_Particle(T_POINT_CLOUD* particle_class)
    {if(particle_class){
        Free_Particle(particle_class->next);particle_class->next=0;particle_class->array_collection->Delete_All_Elements();
#ifdef USE_PTHREADS
        pthread_mutex_lock(&stack_lock);
#endif
        free_pool.Push(particle_class);
#ifdef USE_PTHREADS
        pthread_mutex_unlock(&stack_lock);
#endif
    }}

    T_POINT_CLOUD& Add_Particle(T_POINT_CLOUD& particles,int& index)
    {T_POINT_CLOUD* particles_link=&particles;
    while(particles_link->array_collection->Size()==number_particles_per_cell){ // find the right link (allocate if necessary)
        if(!particles_link->next) particles_link->next=Allocate_Particle();
        particles_link=particles_link->next;}
    index=particles_link->array_collection->Add_Element();
    return *particles_link;}

private:
    void Allocate_New_Batch()
    {PHYSBAM_ASSERT(!template_particles.array_collection->Size()); // make sure our clones are empty
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    LOG::cout<<"Particle Pool: Allocating another "<<allocation_batch_size<<" particles"<<std::endl;
#endif
    CLONE_ARRAY<T_POINT_CLOUD>* batch=new CLONE_ARRAY<T_POINT_CLOUD>(template_particles,allocation_batch_size);
    for(int i=1;i<=allocation_batch_size;i++) (*batch)(i).array_collection->Preallocate(number_particles_per_cell);
    allocated_batches.Append(batch);
#ifdef USE_PTHREADS
    pthread_mutex_lock(&stack_lock);
#endif
    free_pool.Increase_Size(allocation_batch_size);
    for(int i=1;i<=allocation_batch_size;i++) free_pool.Push(&(*batch)(i));
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&stack_lock);
#endif
    }

//#####################################################################
};
}
#endif
