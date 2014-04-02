//#####################################################################
// Copyright 2002-2008, Robert Bridson, Kevin Der, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HASHTABLE
//#####################################################################
#ifndef __HASHTABLE__
#define __HASHTABLE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/PROJECTED_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Math_Tools/Hash.h>
#include <PhysBAM_Tools/Math_Tools/integer_log.h>
#include <PhysBAM_Tools/Utilities/EXCEPTIONS.h>
namespace PhysBAM{

enum HASHTABLE_ENTRY_STATE{ENTRY_FREE,ENTRY_ACTIVE,ENTRY_DELETED};
template<class TK,class T> struct HASHTABLE_ENTRY_TEMPLATE{HASHTABLE_ENTRY_STATE state;TK key;T data;};
template<class TK> struct HASHTABLE_ENTRY_TEMPLATE<TK,void>{HASHTABLE_ENTRY_STATE state;TK key;};

template<class TK,class T> // T = void
class HASHTABLE
{
private:
    struct UNUSABLE{};
    typedef HASHTABLE_ENTRY_TEMPLATE<TK,T> ENTRY; // don't store data if T is void
    typedef typename IF<IS_SAME<T,void>::value,UNUSABLE,T>::TYPE T_UNLESS_VOID;
public:
    typedef TK KEY;
    typedef T ELEMENT;
    typedef ENTRY* iterator; // for stl
    typedef ENTRY value_type; // for stl
    typedef HASHTABLE_ITERATOR<TK,T> ITERATOR;
    typedef HASHTABLE_ITERATOR<TK,const T> CONST_ITERATOR;
    typedef FIELD_PROJECTOR<HASHTABLE_ENTRY_TEMPLATE<TK,T>,HASHTABLE_ENTRY_STATE,&HASHTABLE_ENTRY_TEMPLATE<TK,T>::state> T_FIELD_PROJECTOR;

    ARRAY<ENTRY> table;
    int number_of_entries;
    int next_resize;

private:
    friend class HASHTABLE_ITERATOR<TK,T>;
    friend class HASHTABLE_ITERATOR<TK,const T_UNLESS_VOID>;
public:

    HASHTABLE(const int estimated_max_number_of_entries=5)
    {
        Initialize_New_Table(estimated_max_number_of_entries);
    }

    ~HASHTABLE()
    {}

    void Clean_Memory()
    {Initialize_New_Table(5);} // can't resize to zero since table.m must be a power of two

    int Size() const
    {return number_of_entries;}

    int Max_Size() const
    {return table.m;}

    int Next_Resize() const
    {return next_resize;}

    void Initialize_New_Table(const int estimated_max_number_of_entries_input)
    {next_resize=max(5,estimated_max_number_of_entries_input);
    int estimated_table_entries=(unsigned int)(next_resize*4/3+1); // choose so load is .75
    int number_of_lists=next_power_of_two(estimated_table_entries);
    table.Resize(number_of_lists,false,false); // TODO: only resize if significantly reducing the size
    Remove_All();}

    void Resize_Table(const int estimated_max_number_of_entries_input=0)
    {int estimated_max_number_of_entries=estimated_max_number_of_entries_input;if(!estimated_max_number_of_entries) estimated_max_number_of_entries=3*number_of_entries/2;
    ARRAY<ENTRY> old_table;old_table.Exchange(table);Initialize_New_Table(estimated_max_number_of_entries);
    for(int h=1;h<=old_table.m;h++){ENTRY& entry=old_table(h);if(entry.state==ENTRY_ACTIVE) Insert(entry);}}

private:
    int Next_Index(const int h) const // linear probing
    {return (h&(table.m-1))+1;} // power of two so mod is dropping high order bits

    int Hash_Index(const TK& v) const
    {return (Hash(v)&(table.m-1))+1;} // power of two so mod is dropping high order bits

    bool Contains(const TK& v,const int h) const
    {for(int i=h;;i=Next_Index(i))
        if(table(i).state==ENTRY_FREE) return false;
        else if(table(i).state==ENTRY_ACTIVE && table(i).key==v) return true;}

    void Insert(const HASHTABLE_ENTRY_TEMPLATE<TK,void>& entry)
    {Insert(entry.key);}

    template<class T2> typename DISABLE_IF<IS_SAME<T2,void>::value>::TYPE
    Insert(const HASHTABLE_ENTRY_TEMPLATE<TK,T2>& entry)
    {Insert(entry.key,entry.data);}

public:

    void Insert(const TK& v) // assumes no entry with v exists
    {STATIC_ASSERT((IS_SAME<T,void>::value));
    if(number_of_entries>next_resize) Resize_Table();
    number_of_entries++;
    int h=Hash_Index(v);assert(!Contains(v,h));
    for(;table(h).state!=ENTRY_FREE;h=Next_Index(h));
    ENTRY& entry=table(h);entry.key=v;entry.state=ENTRY_ACTIVE;}

    T_UNLESS_VOID& Insert(const TK& v,const T_UNLESS_VOID& value) // assumes no entry with v exists
    {STATIC_ASSERT((NOT<IS_SAME<T,void>::value>::value));
    if(number_of_entries>next_resize) Resize_Table();
    number_of_entries++;
    int h=Hash_Index(v);assert(!Contains(v,h));
    for(;table(h).state!=ENTRY_FREE;h=Next_Index(h));
    ENTRY& entry=table(h);entry.key=v;entry.state=ENTRY_ACTIVE;entry.data=value;return entry.data;}

    T_UNLESS_VOID& Get_Or_Insert(const TK& v,const T_UNLESS_VOID& default_value=T_UNLESS_VOID()) // inserts the default if key not found
    {int h=Hash_Index(v);
    for(;table(h).state!=ENTRY_FREE;h=Next_Index(h))
        if(table(h).state==ENTRY_ACTIVE && table(h).key==v) return table(h).data;
    if(number_of_entries>next_resize){return Insert(v, default_value);}
    number_of_entries++;
    ENTRY& entry=table(h);entry.key=v;entry.state=ENTRY_ACTIVE;entry.data=default_value;return entry.data;}

    T* Get_Pointer(const TK& v) // returns NULL if key not found
    {for(int h=Hash_Index(v);table(h).state!=ENTRY_FREE;h=Next_Index(h))
        if(table(h).state==ENTRY_ACTIVE && table(h).key==v) return &table(h).data;
    return 0;}

    const T* Get_Pointer(const TK& v) const // returns NULL if key not found
    {for(int h=Hash_Index(v);table(h).state!=ENTRY_FREE;h=Next_Index(h))
        if(table(h).state==ENTRY_ACTIVE && table(h).key==v) return &table(h).data;
    return 0;}

    T_UNLESS_VOID& Get(const TK& v) // fails if key not found
    {if(T_UNLESS_VOID* data=Get_Pointer(v)) return *data;throw KEY_ERROR("HASHTABLE::Get");}

    const T_UNLESS_VOID& Get(const TK& v) const // fails if key not found
    {if(const T_UNLESS_VOID* data=Get_Pointer(v)) return *data;throw KEY_ERROR("HASHTABLE::Get");}

    T Get_Default(const TK& v,const T_UNLESS_VOID& default_value=T_UNLESS_VOID()) const // returns default_value if key not found
    {if(const T* data=Get_Pointer(v)) return *data;return default_value;}

    bool Contains(const TK& v) const
    {return Contains(v,Hash_Index(v));}

    bool Get(const TK& v,T_UNLESS_VOID& value) const
    {for(int h=Hash_Index(v);table(h).state!=ENTRY_FREE;h=Next_Index(h))
         if(table(h).state==ENTRY_ACTIVE && table(h).key==v){value=table(h).data;return true;}
    return false;}

    bool Set(const TK& v) // insert entry if doesn't already exists, returns whether it added a new entry
    {STATIC_ASSERT((IS_SAME<T,void>::value));
    if(number_of_entries>next_resize) Resize_Table();
    int h=Hash_Index(v);
    for(;table(h).state!=ENTRY_FREE;h=Next_Index(h))
        if(table(h).state==ENTRY_ACTIVE && table(h).key==v) return false;
    number_of_entries++;
    ENTRY& entry=table(h);entry.key=v;entry.state=ENTRY_ACTIVE;return true;}

    bool Set(const TK& v,const T_UNLESS_VOID& value) // if v doesn't exist insert value, else sets its value, returns whether it added a new entry
    {STATIC_ASSERT((NOT<IS_SAME<T,void>::value>::value));
    if(number_of_entries>next_resize) Resize_Table(); // if over load average, have to grow (must do this before computing hash index)
    int h=Hash_Index(v);
    for(;table(h).state!=ENTRY_FREE;h=Next_Index(h))
        if(table(h).state==ENTRY_ACTIVE && table(h).key==v){table(h).data=value;return false;}
    number_of_entries++;
    ENTRY& entry=table(h);entry.key=v;entry.state=ENTRY_ACTIVE;entry.data=value;return true;}

    template<class T_ARRAY>
    void Set_All(const T_ARRAY& array)
    {STATIC_ASSERT(IS_SAME<typename T_ARRAY::ELEMENT,TK>::value);
    for(typename T_ARRAY::INDEX i(1);i<=array.Size();i++) Set(array(i));}

    bool Delete_If_Present(const TK& v)
    {for(int h=Hash_Index(v);table(h).state!=ENTRY_FREE;h=Next_Index(h)) // reduce as still are using entries for deletions
        if(table(h).state==ENTRY_ACTIVE && table(h).key==v){table(h).state=ENTRY_DELETED;number_of_entries--;next_resize--;return true;}
    return false;}

    void Delete(const TK& v)
    {if(!Delete_If_Present(v)) throw KEY_ERROR("HASHTABLE::Delete");}

    void Remove_All()
    {PROJECTED_ARRAY<ARRAY<HASHTABLE_ENTRY_TEMPLATE<TK,T>,int>,T_FIELD_PROJECTOR> projected_array=table.template Project<HASHTABLE_ENTRY_STATE,&HASHTABLE_ENTRY_TEMPLATE<TK,T>::state>();
        ARRAYS_COMPUTATIONS::Fill(projected_array,ENTRY_FREE);number_of_entries=0;}
    //{table.template Project<HASHTABLE_ENTRY_STATE,&HASHTABLE_ENTRY_TEMPLATE<TK,T>::state>().Fill(ENTRY_FREE);number_of_entries=0;}

    void Exchange(const TK& x,const TK& y) // Exchange values at entries x and y; valid if x or y (or both) are not present; efficient for array values
    {bool a=Contains(x),b=Contains(y);if(a || b){exchange(Get_Or_Insert(x),Get_Or_Insert(y));if(!a || !b){Delete(a?x:y);}}}

    static void Exchange_Hashtables(HASHTABLE& hash1,HASHTABLE& hash2)
    {hash1.table.Exchange(hash2.table);exchange(hash1.number_of_entries,hash2.number_of_entries);exchange(hash1.next_resize,hash2.next_resize);}

    void Print_Table(std::ostream& output) const
    {output<<"Entry Count: "<<number_of_entries<<" Elements to resize at: "<<next_resize<<std::endl;
    for(int h=1;h<=table.m;h++)
         output<<h<<":"<<(table(h).state==ENTRY_ACTIVE?"ACTIVE":(table(h).state==ENTRY_DELETED?"DELETED":"FREE"))<<" key="<<table(h).key<<" value="<<table(h).data<<std::endl;}

    void Apply_Function_To_All_Entries(void (*function)(TK&,T_UNLESS_VOID&))
    {for(int h=1;h<=table.m;h++) if(table(h).state==ENTRY_ACTIVE) function(table(h).key,table(h).data);}

    void Delete_Pointers_Stored_In_Table() // of course, only valid if pointers are stored in table
    {for(int h=1;h<=table.m;h++) if(table(h).state==ENTRY_ACTIVE){delete table(h).data;table(h).data=0;}}

    void Reset_List_Arrays_Stored_In_Table() // of course, only works if pointers to ARRAY are stored in table
    {for(int h=1;h<=table.m;h++) if(table(h).state==ENTRY_ACTIVE){table(h).data->Remove_All();}}

    void Get_Keys(ARRAY<TK>& keys) const
    {keys.Remove_All();keys.Preallocate(Size());for(int h=1;h<=table.m;h++) if(table(h).state==ENTRY_ACTIVE) keys.Append(table(h).key);}

    template<class T_ARRAY1,class T_ARRAY2>
    void Get_Complementary_Keys(const T_ARRAY1& keys_universe,T_ARRAY2& keys_complementary) const
    {STATIC_ASSERT(IS_SAME<typename T_ARRAY1::ELEMENT,TK>::value && IS_SAME<typename T_ARRAY2::ELEMENT,TK>::value);
    keys_complementary.Remove_All();keys_complementary.Preallocate(keys_universe.Size()-Size());
    for(typename T_ARRAY1::INDEX i(1);i<=keys_universe.Size();i++) if(!Contains(keys_universe(i))) keys_complementary.Append(keys_universe(i));}

    void Get_Data(ARRAY<T_UNLESS_VOID>& data) const
    {data.Remove_All();data.Preallocate(Size());for(int h=1;h<=table.m;h++) if(table(h).state==ENTRY_ACTIVE) data.Append(table(h).data);}

    ENTRY* begin()
    {return table.Get_Array_Pointer();}

    const ENTRY* begin() const
    {return table.Get_Array_Pointer();}

    ENTRY* end()
    {return table.Get_Array_Pointer()+table.Size();}

    const ENTRY* end() const
    {return table.Get_Array_Pointer()+table.Size();}

//#####################################################################
};
}
#endif
