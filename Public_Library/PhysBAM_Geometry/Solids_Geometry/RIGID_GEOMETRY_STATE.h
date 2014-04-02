//#####################################################################
// Copyright 2003-2006, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GEOMETRY_STATE
//#####################################################################
#ifndef __RIGID_GEOMETRY_STATE__
#define __RIGID_GEOMETRY_STATE__

#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
namespace PhysBAM{

template<class TV>
class RIGID_GEOMETRY_STATE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
public:
    T time;
    FRAME<TV> frame;
    TWIST<TV> twist;

    RIGID_GEOMETRY_STATE()
        :time(0)
    {}

    RIGID_GEOMETRY_STATE(const FRAME<TV>& frame_input)
       :time(0),frame(frame_input)
    {}

    RIGID_GEOMETRY_STATE(const FRAME<TV>& frame_input,const TWIST<TV>& twist_input)
       :time(0),frame(frame_input),twist(twist_input)
    {}

    template<class TV2> explicit RIGID_GEOMETRY_STATE(const RIGID_GEOMETRY_STATE<TV2>& state_input)
        :time((T)state_input.time),frame(state_input.frame),twist(state_input.twist)
    {}

    static void Compute_Velocity_Between_States(const RIGID_GEOMETRY_STATE& state1,const RIGID_GEOMETRY_STATE& state2,RIGID_GEOMETRY_STATE& result_state)
    {T one_over_dt=1/(state2.time-state1.time);
    result_state.twist.linear=one_over_dt*(state2.frame.t-state1.frame.t);
    result_state.twist.angular=one_over_dt*(state2.frame.r*state1.frame.r.Inverse()).Rotation_Vector();}

    TV Object_Space_Point(const TV& world_space_point) const
    {return frame.Inverse_Times(world_space_point);}

    TV Object_Space_Vector(const TV& world_space_vector) const
    {return frame.r.Inverse_Rotate(world_space_vector);}

    TV World_Space_Point(const TV& object_space_point) const
    {return frame*object_space_point;}

    TV World_Space_Vector(const TV& object_space_vector) const
    {return frame.r.Rotate(object_space_vector);}

    TV Pointwise_Object_Velocity(const TV& X) const
    {return twist.linear+TV::Cross_Product(twist.angular,X-frame.t);}

    template<class RW> static void Read_Frame(std::istream& input,FRAME<VECTOR<T,1> >& f)
    {Read_Binary<RW>(input,f);}

    template<class RW> static void Read_Frame(std::istream& input,FRAME<VECTOR<T,2> >& f)
    {VECTOR<T,1> angle;Read_Binary<RW>(input,f.t,angle);f.r=ROTATION<VECTOR<T,2> >::From_Rotation_Vector(angle);}

    template<class RW> static void Read_Frame(std::istream& input,FRAME<VECTOR<T,3> >& f)
    {Read_Binary<RW>(input,f);}

    template<class RW> static void Write_Frame(std::ostream& output,const FRAME<VECTOR<T,1> >& f)
    {Write_Binary<RW>(output,f);}

    template<class RW> static void Write_Frame(std::ostream& output,const FRAME<VECTOR<T,2> >& f)
    {Write_Binary<RW>(output,f.t,f.r.Angle());}

    template<class RW> static void Write_Frame(std::ostream& output,const FRAME<VECTOR<T,3> >& f)
    {Write_Binary<RW>(output,f);}

//#####################################################################
};
}
#include <PhysBAM_Geometry/Read_Write/Solids_Geometry/READ_WRITE_RIGID_GEOMETRY_STATE.h>
#endif
