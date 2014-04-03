//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS
//#####################################################################
#ifndef __COMBINED_COLLISIONS__
#define __COMBINED_COLLISIONS__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_BODY_ID.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

static const int COMBINED_COLLISIONS_TYPE_RIGID=0x20000000;
static const int COMBINED_COLLISIONS_TYPE_DEFORMABLE=0x40000000;

inline bool Combined_Body_Id_Is_Rigid(COMBINED_BODY_ID id)
{
    return Value(id)&COMBINED_COLLISIONS_TYPE_RIGID;
}

inline int Combined_Body_Id_To_Rigid_Body(COMBINED_BODY_ID id)
{
    assert(Value(id)&COMBINED_COLLISIONS_TYPE_RIGID);
    return Value(id)&~COMBINED_COLLISIONS_TYPE_RIGID;
}

inline int Combined_Body_Id_To_Deformable_Particle(COMBINED_BODY_ID id)
{
    assert(Value(id)&COMBINED_COLLISIONS_TYPE_DEFORMABLE);
    return Value(id)&~COMBINED_COLLISIONS_TYPE_DEFORMABLE;
}

inline COMBINED_BODY_ID Rigid_Body_To_Combined_Body_Id(int p)
{
    return COMBINED_BODY_ID(p|COMBINED_COLLISIONS_TYPE_RIGID);
}

inline COMBINED_BODY_ID Deformable_Particle_To_Combined_Body_Id(int p)
{
    return COMBINED_BODY_ID(p|COMBINED_COLLISIONS_TYPE_DEFORMABLE);
}

template<class TV>
class COMBINED_COLLISIONS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
public:
//    RIGID_DEFORMABLE_COLLISIONS<TV>& rigid_deformable_collisions;

    struct IMPULSE
    {
        IMPULSE();
        virtual ~IMPULSE();

        virtual void Resize()=0;
        virtual void Apply(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time)=0;
        virtual void Clear(const ARRAY<COMBINED_BODY_ID>& list)=0;
        virtual void Clear()=0;
        virtual void Setup_Discover_State(T dt,T time)=0;
        virtual void Setup_Impulse_State(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time)=0;
        virtual void Finish_Step(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time)=0;
        virtual T Kinetic_Energy(const ARRAY<COMBINED_BODY_ID>& list) const=0;
    };

    struct COLLIDER:public NONCOPYABLE
    {
        COLLIDER();
        virtual ~COLLIDER();

        virtual void Discover(const T dt,const T time)=0;
        virtual int Count() const=0;
        virtual void Precompute(int e)=0;
        virtual void Affected_Bodies(int e,HASHTABLE<COMBINED_BODY_ID>& hash) const=0;
        virtual bool Accumulate_Impulse(int e,IMPULSE& impulse,T scale,T dt,T time) const=0;
        virtual T Diagnose_Impulse(int e,const IMPULSE& impulse,T dt,T time) const=0;
        virtual void Moved(int e)=0;
        virtual T Kinetic_Energy_Gradient(int e,const IMPULSE& impulse,T scale,T dt,T time) const=0;
    };

    struct COLLIDER_INFO
    {
        ARRAY<T> scale;
    };

    IMPULSE* accumulator;
    ARRAY<COLLIDER*> colliders;
    ARRAY<COLLIDER_INFO> collider_info;
    HASHTABLE<COMBINED_BODY_ID> referenced_particles;
    ARRAY<COMBINED_BODY_ID> referenced_list;
    bool test_system;

    COMBINED_COLLISIONS(IMPULSE* accumulator_input);
    ~COMBINED_COLLISIONS();

    void Add_Collider(COLLIDER* collider);
    void Apply_Combined_Impulse(T dt,T time);
    void Apply_Combined_Projections(T dt,T time);
    void Update_Scales_With_Diagnosis(T dt,T time);
    void Refresh_Referenced_List();
    void Compute_Impulse_Using_Scale(T dt,T time);
    void Apply_Computed_Impulses(T dt,T time);
    void Notify_Moved();
    void Solve(T dt,T time,int iterations);
    void Iterate_Projections(T dt,T time,int iterations);
    void Recompute_Active_Set(T dt,T time);

    void Test_System();
    void Find_Scales_With_CG(T dt,T time);
//#####################################################################
};
}
#endif
