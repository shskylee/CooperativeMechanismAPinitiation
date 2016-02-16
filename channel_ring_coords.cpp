
//  Created by ShanshanLi on 7/21/15.
//
//

#include "channel_ring_coords.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
#include <iostream>

channel_ring_coords::channel_ring_coords(){
    D_clus=0.190;
    D_ring=0.031;
    AD_chl=0.1653469;
    //Nclus=158;
    
    // test algorithm
    Nclus=2;
    // set cluster number as 2
    
    Nring=2;
    Nrchl=38;
}

channel_ring_coords::channel_ring_coords(const int u, const int v , const int w, const double x, const double y, const double z){
    Nclus=u;
    Nring=v;
    Nrchl=w;
    D_clus=x;
    D_ring=y;
    AD_chl=z;
}

boost::numeric::ublas::matrix<double> channel_ring_coords::ring_coordinates()
{

    int Nchl=Nclus*Nring*Nrchl;
    int channel;
    boost::numeric::ublas::matrix<double> ring_coordinates(Nchl,3);
    for (int i=0; i<Nclus; i++) {
        for (int j=0; j<Nring; j++) {
            for (int k=0; k<Nrchl; k++) {
                channel=i*Nring*Nrchl+j*Nrchl+k;
                ring_coordinates(channel,0)=i*D_clus+j*D_ring;
                ring_coordinates(channel,1)=cos(AD_chl*k);
                ring_coordinates(channel,2)=sin(AD_chl*k);
            }
        }
    }
    return ring_coordinates;
}

