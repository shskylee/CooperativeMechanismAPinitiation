//
//  dependency_graph.cpp
//  
//
//  Created by ShanshanLi on 7/21/15.
//
//

#include "channel_coordinate.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

boost::numeric::ublas::matrix<double> channel_coordinate::coordinates()
{
    // distances between Nchl channels
    boost::mt19937 rng;
    boost::uniform_real<> uni_dist(-1.0*Radius,1.0*Radius);
    boost::variate_generator<boost::mt19937&,boost::uniform_real<double> > gen_random(rng,uni_dist);
    boost::numeric::ublas::matrix<double> coordinate(Nchl,3);
    for (int i=0; i<coordinate.size1(); i++) {
        for (int j=0; j<coordinate.size2(); j++) {
           coordinate(i,j)=gen_random();
        }
    }
    return coordinate;
}

channel_coordinate::channel_coordinate(const int x, const double y){
    Nchl=x;
    Radius=y;
}