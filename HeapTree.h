//
//  heapTree.h
//  
//
//  Created by ShanshanLi on 7/27/15.
//
//

#ifndef _HeapTree_h
#define _HeapTree_h


#include <assert.h>  // for error-checking purposes

template <class Elem>
class HeapTree
{
public:
    HeapTree(int MaxSize=100000);
    HeapTree(const HeapTree<Elem> &OtherTree);
    HeapTree(Elem *Array, int ElemNum, int MaxSize);
    Elem *Sort(void); // Built-in HeapSort Algorithm
    ~HeapTree(void);
    
    bool Add(const Elem &Item); //Add the Item to Heap
    Elem Remove(void);         //Remove and return Item from Heap
    bool Update(int Node,const Elem &value);  // upate value of position to new value
    bool Update_Aux(int Node); // update aux func
    bool Root_Node(Elem &root_data,int &root_index);
    
    inline int GetSize(void);  //Returns the number of nodes in the heap
    
private:
    Elem  *Data; // actual data array
    int   *Data_Index; // index of data
    int   *Index_Node;    // ith element is a pointer to the position on the tree of ith data point
    int    CurrentNum; //Current number of elements
    const int MAX_SIZE; // the maximum number of elements
    
    void ShiftUp (int Node);  //Shift node up into place
    void ShiftDown (int Node); //Shift node down into place
    
    inline int ParentOf(int Node);   //Returns Parent location
    inline int LeftChildOf(int Node);   //Returns Left Child location
};

//HeapTree constructor function
template <class Elem>
HeapTree<Elem>::HeapTree(int MaxSize):MAX_SIZE(MaxSize)
{
    Data = new Elem[MAX_SIZE];
    Data_Index= new int[MAX_SIZE];
    Index_Node= new int[MAX_SIZE];
    CurrentNum=0;
}

//HeapTree copy constructor function
template<class Elem>
HeapTree<Elem>::HeapTree(const HeapTree<Elem> &OtherTree):MAX_SIZE(OtherTree.MAX_SIZE)
{
    Data = new Elem[MAX_SIZE];
    Data_Index = new int[MAX_SIZE];
    Index_Node= new int[MAX_SIZE];
    
    CurrentNum = OtherTree.CurrentNum;
    
    //copy the array
    for (int i=0; i<OtherTree.CurrentNum; ++i) {
        Data[i]=OtherTree.Data[i];
        Data_Index[i] = OtherTree.Data_Index[i];
        Index_Node[i]= OtherTree.Index_Node[i];
    }
}

//HeapTree array constructor
template<class Elem>
HeapTree<Elem>::HeapTree(Elem *Array, int ElemNum, int MaxSize): MAX_SIZE(MaxSize)
{
    Data = new Elem[MAX_SIZE];
    Data_Index = new int[MAX_SIZE];
    Index_Node= new int[MAX_SIZE];
    CurrentNum = ElemNum;
    
    // this organizes the array into a proper HeapTree
    for (int i=0; i<ElemNum; ++i) {
        Data[i]=Array[i];
        Data_Index[i]=i;
        Index_Node[i]=i;
    }
    // this organizes the array into a proper HeapTree
    for (int i=ParentOf(CurrentNum-1); i>=0; --i) {
        ShiftDown(i);
    }
}
// Built-in Heap Sort algorithm
template <class Elem>
Elem *HeapTree<Elem>::Sort(void)
{
    // This is the array that will be returned
    Elem *NewArray = new Elem[CurrentNum];
    
    // The algorithm works back to front, with the sorted
    // elements being stored in NewArray
    for (int ElemNum = 0; ElemNum <CurrentNum; ++ElemNum)
    {
        // Since the Remove() function alters CurrentNum by subtracting 1
        // from it each time, we must use a seperate variable to
        // index NewArray.
        NewArray[ElemNum] = Remove();
    }
    return NewArray;
}
// return the data and index of the first node
template <class Elem>
bool HeapTree<Elem>::Root_Node(Elem &root_data, int &root_index)
{
    root_data=Data[0];
    root_index=Data_Index[0];
    return true;
}

// HeapTree desctructor function
template <class Elem>
HeapTree<Elem>::~HeapTree(void)
{
    if (Data) {
        delete Data;
    }
    if (Data_Index) {
        delete Data_Index;
    }
    if (Index_Node) {
        delete Index_Node;
    }
}

// Add() function
template <class Elem>
bool HeapTree<Elem>::Add(const Elem &Item)
{
    if (CurrentNum >= MAX_SIZE) {
        return false;    // if we have reached our maximum capacity
    }
    Data[CurrentNum] = Item;
    Data_Index[CurrentNum]=CurrentNum;
    Index_Node[CurrentNum]=CurrentNum;
    ShiftUp(CurrentNum++);
    return true;
}
// remove() function
template <class Elem>
Elem HeapTree<Elem>::Remove(void)
{
    assert(CurrentNum>0);
    
    Elem Temp=Data[0];
    Data[0]=Data[--CurrentNum]; //replace with the last element
    Data_Index[0]=Data_Index[CurrentNum];
    Index_Node[Data_Index[0]]=0;
    ShiftDown(0);
    
    return Temp;
    
}
template <class Elem>
bool HeapTree<Elem>::Update(int Node,const Elem &value)
{
    int Current=Node;
    if (Current>=CurrentNum) {
        return false;  // such position exceeds current capacity
    }
    Data[Current]=value;
    Update_Aux(Current);
    return true;
}
template <class Elem>
bool HeapTree<Elem>::Update_Aux(int Node)
{
    int Current=Node,
        Parent  = ParentOf(Current),
        Child   = LeftChildOf(Current);
    
    if (Current>=CurrentNum) {
        return false;
    }
    
    if (Child < (CurrentNum - 1))
        if (Data[Child] > Data[Child+1])  // Set Child to smallest Child node
            ++Child;
    
    if (Data[Current]<Data[Parent])
        ShiftUp(Current);
    else if(Data[Current]>Data[Child])
        ShiftDown(Current);
    else
        return true;
}

// GetSize() function
template<class Elem>
inline int HeapTree<Elem>::GetSize(void)
{
    return CurrentNum;
}
// ShiftUp() function
template <class Elem>
void HeapTree<Elem>::ShiftUp(int Node)
{
    int  Current = Node,
    Parent  = ParentOf(Current);
    Elem Item    = Data[Current];
    int  Index  = Data_Index[Current];  // index of reaction at the current position
    
    while (Current > 0)  // While Current is not the RootNode
    {
        if (Data[Parent] > Item)
        {
            Data[Current] = Data[Parent];
            Data_Index[Current]=Data_Index[Parent];
            Index_Node[Data_Index[Current]]=Current;
            Current = Parent;
            Parent = ParentOf(Current);
        }
        else
            break;
    }
    Data[Current] = Item;
    Data_Index[Current]=Index;
    Index_Node[Index]=Current;
}

// ShiftDown() function
template <class Elem>
void HeapTree<Elem>::ShiftDown(int Node)
{
    int Current = Node,
    Child   = LeftChildOf(Current);
    Elem Item   = Data[Current];    // Used to compare values
    int  Index  = Data_Index[Current];  // index of reaction at the current node
   
    
    while (Child < CurrentNum)
    {
        if (Child < (CurrentNum - 1))
            if (Data[Child] > Data[Child+1])  // Set Child to smallest Child node
                ++Child;
        
        if (Item > Data[Child])
        {    // Switch the Current node and the Child node
            Data[Current] = Data[Child];
            Data_Index[Current]=Data_Index[Child];
            Index_Node[Data_Index[Current]]=Current;
            Current       = Child;
            Child         = LeftChildOf(Current);
        }
        else
            break;
    }
    Data[Current] = Item;
    Data_Index[Current]=Index;
    Index_Node[Index]=Current;
}
// ParentOf() function
template <class Elem>
inline int HeapTree<Elem>::ParentOf(int Node)
{
    assert(Node > 0);
    // This uses the fact that decimals are truncated during
    // the division of integers. Thus, (12 - 1) / 2 == 5
    return (Node - 1) / 2;
}

// LeftChildOf() function
template <class Elem>
inline int HeapTree<Elem>::LeftChildOf(int Node)
{
    return (Node * 2) + 1;
}

#endif
