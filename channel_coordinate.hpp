//
//  dependency_graph.hpp
//  
//
//  Created by ShanshanLi on 7/21/15.
//
//

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>

class channel_coordinate
{
private:
    double Radius; // Radius of 3-d space of channel distribution
    int    Nchl;   // number of channels
public:
    //constructors
    channel_coordinate(int,double);
    boost::numeric::ublas::matrix<double>  coordinates();
   
};