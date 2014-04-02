//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_HASHTABLE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_HASHTABLE__
#define __READ_WRITE_HASHTABLE__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class TK,class T>
class Read_Write<HASHTABLE<TK,T>,RW>
{
    struct UNUSABLE{};
public:
    static void Read(typename IF<IS_SAME<T,void>::value,std::istream&,UNUSABLE>::TYPE input,HASHTABLE<TK,T>& object) // void version
    {int entries;Read_Binary<RW>(input,entries);object.Initialize_New_Table(entries);
    for(int i=1;i<=entries;i++){TK key;Read_Binary<RW>(input,key);object.Insert(key);}}

    static void Read(typename IF<IS_SAME<T,void>::value,UNUSABLE,std::istream&>::TYPE input,HASHTABLE<TK,T>& object) // non-void version
    {int entries;Read_Binary<RW>(input,entries);object.Initialize_New_Table(entries);
    for(int i=1;i<=entries;i++){TK key;T value;Read_Binary<RW>(input,key,value);object.Insert(key,value);}}

    static void Write(typename IF<IS_SAME<T,void>::value,std::ostream&,UNUSABLE>::TYPE output,const HASHTABLE<TK,T>& object) // void version
    {Write_Binary<RW>(output,object.number_of_entries);
    for(int h=1;h<=object.table.m;h++) if(object.table(h).state==ENTRY_ACTIVE) Write_Binary<RW>(output,object.table(h).key);}

    static void Write(typename IF<IS_SAME<T,void>::value,UNUSABLE,std::ostream&>::TYPE output,const HASHTABLE<TK,T>& object) // non-void version
    {Write_Binary<RW>(output,object.number_of_entries);
    for(int h=1;h<=object.table.m;h++) if(object.table(h).state==ENTRY_ACTIVE) Write_Binary<RW>(output,object.table(h).key,object.table(h).data);}
};

template<class K,class T>
std::ostream& operator<<(std::ostream& output,const HASHTABLE<K,T>& hashtable)
{
    output<<"(";
    bool first=true;
    for(typename HASHTABLE<K,const T>::ITERATOR iterator(hashtable);iterator.Valid();iterator.Next()){
        if(!first) output<<" ";first=false;
        output<<iterator.Key()<<":"<<iterator.Data();}
    output<<")";
    return output;
}

template<class K>
std::ostream& operator<<(std::ostream& output,const HASHTABLE<K,void>& hashtable)
{
    output<<"(";
    bool first=true;
    for(typename HASHTABLE<K,void>::ITERATOR iterator(hashtable);iterator.Valid();iterator.Next()){
        if(!first) output<<" ";first=false;
        output<<iterator.Key();}
    output<<")";
    return output;
}
}

#endif
#endif
