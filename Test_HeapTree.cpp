//
//  Test_queue.cpp
//  
//
//  Created by ShanshanLi on 7/27/15.
//
//

#include "HeapTree.h"
#include <iostream>
#include <fstream>
#include <iomanip>


int main(int argc, char **argv)
{
    
    double Array[10]={4.2,3.7,1000,5.5,9.1,2.0,8.9,1.3,7.3,10.1};

    HeapTree<double> HeapT(Array, 10, 20);
    
  
    
    

        std::cout<<HeapT.Sort()<<"\n";
    
    return 0;
}
