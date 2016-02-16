//
//  channel_ring_coords.hpp
//  we assume the ring distribution
//  number of clusters: 158
//  number of channels per cluster: 76
//  number of rings per cluster: 2
//  number of channels per ring: 38
//  angular distance between channels in one cluster :2*pi/38
//  distance between clusters: 190 nm 0.19 micrometers

// the diameter of axon initial segment: 1 micrometer
// the length of axon initial segment: 30 micrometers

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


class channel_ring_coords
{
private:
    double D_clus;  // Distance beteen clusters
    double D_ring;  // Distance between rings within one cluster
    double AD_chl;  // angular distance between channels within one ring
    int    Nclus;   // number of clusters
    int    Nring;   // number of rings per clusters
    int    Nrchl;   // number of channels per ring
public:
    //constructors with default parameters
    channel_ring_coords();
    //constructors with given parameters
    channel_ring_coords(int,int,int,double,double,double);
    boost::numeric::ublas::matrix<double>  ring_coordinates();
};