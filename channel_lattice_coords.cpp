//
//  dependency_graph.cpp
//  
//
//  Created by ShanshanLi on 7/21/15.
//
//

#include "channel_lattice_coords.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
#include <iostream>

channel_lattice_coords::channel_lattice_coords(){
    D_clus=2.0;
    D_chl=0.6;
    Nclus=20;
    Npclus=100;
}

channel_lattice_coords::channel_lattice_coords(const int x, const int x0 ,const double y, const double z){
    Nclus=x;
    Npclus=x0;
    D_clus=y;
    D_chl=z;
}

boost::numeric::ublas::matrix<double> channel_lattice_coords::lattice_coordinates()
{

    int Nchl=Nclus*Npclus;
    int m=sqrt(Npclus);
    int channel;
    boost::numeric::ublas::matrix<double> lattice_coordinates(Nchl,3);
    for (int i=0; i<Nclus; i++) {
        for (int j=0; j<m; j++) {
            for (int k=0; k<m; k++) {
                channel=i*Npclus+j*m+k;
                lattice_coordinates(channel,0)=i*D_clus;
                lattice_coordinates(channel,1)=j*D_chl;
                lattice_coordinates(channel,2)=k*D_chl;
            }
        }
    }
    return lattice_coordinates;
}

