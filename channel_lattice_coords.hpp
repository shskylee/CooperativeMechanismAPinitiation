//
//  channel_lattice_coords.hpp
//  we assume the lattice distribution
//  number of cluster 20
//  number of channels per cluster 100
//  distance between clusters 2
//  distance between channels within one cluster 0.6


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


class channel_lattice_coords
{
public:
    double D_clus; // Distance beteen channels
    double D_chl;  // closest distance between channels within one cluster
    int    Nclus;   // number of clusters
    int    Npclus;  // number of clsuter
public:
    //constructors with default parameters
    channel_lattice_coords();
    //constructors with given parameters
    channel_lattice_coords(int,int,double,double);
    boost::numeric::ublas::matrix<double>  lattice_coordinates();
};