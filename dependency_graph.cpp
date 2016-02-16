//
//  dependency_graph.cpp
//  
//
//  Created by ShanshanLi on 7/24/15.
//
//

#include "dependency_graph.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <fstream>

dependency_graph::dependency_graph(int u,int v, double x)
{
    Nchl=u;
    Nb=v;
    diffusion_thres=x;
    
}

boost::numeric::ublas::matrix<double> dependency_graph::channel_neighbors()
{
 
    boost::numeric::ublas::matrix<double>  neighbors(Nb*Nchl,2); // neighbors for all channels
    std::string line;
    std::ifstream chl_neighbor;
    chl_neighbor.open("channel_neighbors.txt");
    int i=0, j;
    double a;
    while (getline(chl_neighbor, line) && !line.empty()) {
        std::istringstream iss(line);
        j = 0;
        while (iss >> a) {
            neighbors(i,j) = a;
            j++;
        }
        i++;
    }

    chl_neighbor.close();
    return neighbors;
    
}

boost::numeric::ublas::matrix<int>  dependency_graph::Depend_graph()
{
    // dependency graph with chosen reaction diagram without considering the local voltage
    //number of reactions to one channel, initially set to 9
    //int Nres=9;
    int Nres=5;
    int i,j,k;
    int neighbor;
    boost::numeric::ublas::matrix<double>  neighbors(Nb*Nchl,2); // neighbors for all channels
    boost::numeric::ublas::zero_matrix<int> m(Nres, Nres);
    boost::numeric::ublas::matrix<int>      Demat=m;
    
    //chosen reaction diagram in order : c->m (0),c->i (1),m->c(2),i->c(3),m->i(4),i->m(5),m->o(6),o->m(7),o->i(8)
    /*for (i=0; i<Nres; i++) {
        Demat(i,i)=1;
    }
    Demat(0,1)=1;Demat(0,2)=1;Demat(0,4)=1;Demat(0,6)=1;
    Demat(1,0)=1;Demat(1,3)=1;Demat(1,5)=1;
    Demat(2,0)=1;Demat(2,1)=1;Demat(2,4)=1;Demat(2,6)=1;
    Demat(3,0)=1;Demat(3,1)=1;Demat(3,5)=1;
    Demat(4,2)=1;Demat(4,3)=1;Demat(4,5)=1;Demat(4,6)=1;
    Demat(5,2)=1;Demat(5,3)=1;Demat(5,4)=1;Demat(5,6)=1;
    Demat(6,2)=1;Demat(6,4)=1;Demat(6,7)=1;Demat(6,8)=1;
    Demat(7,2)=1;Demat(7,4)=1;Demat(7,6)=1;Demat(7,8)=1;
    Demat(8,3)=1;Demat(8,5)=1;Demat(8,7)=1;*/
    
    //chosen reaction diagram in order: c->o (0), o->c (1), i->c (2), c->i(3), o->i (4)
    
    for (i=0; i<Nres; i++){
        Demat(i,i)=1;
    }
    
    Demat(0,1)=1;Demat(0,3)=1;Demat(0,4)=1;
    Demat(1,0)=1;Demat(1,3)=1;Demat(1,4)=1;
    Demat(2,0)=1;Demat(2,3)=1;
    Demat(3,0)=1;Demat(3,2)=1;
    Demat(4,1)=1;Demat(4,2)=1;
  

    neighbors=channel_neighbors();
    boost::numeric::ublas::zero_matrix<int> m1(Nchl*Nres, Nchl*Nres);
    boost::numeric::ublas::matrix<int>  Gmat=m1;
    for (i=0; i<Nchl; i++) {
        boost::numeric::ublas::project(Gmat,boost::numeric::ublas::range(i*Nres,i*Nres+Nres),boost::numeric::ublas::range(i*Nres,i*Nres+Nres))=Demat;
    }
    
    i=0;
    while (i<Nchl) {
         for (j=0; j<Nb; j++) {
             if (neighbors(i*Nb+j,1)<=diffusion_thres) {
                 neighbor=neighbors(i*Nb+j,0);
                 for (k=0; k<Nres-1; k++) {
                     Gmat(i*Nres,neighbor*Nres+k)=1;
                     Gmat(i*Nres+1,neighbor*Nres+k)=1;
                     Gmat(i*Nres+4,neighbor*Nres+k)=1;
                 }
                /* for (k=0; k<Nres; k++) {
                     Gmat(i*Nres+6,neighbor*Nres+k)=1;
                     Gmat(i*Nres+7,neighbor*Nres+k)=1;
                     Gmat(i*Nres+8,neighbor*Nres+k)=1;
                     
                 }*/
             }
         }
        i++;
    }
    
    return Gmat;
   
}
