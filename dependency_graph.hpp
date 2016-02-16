//
//  dependency_graph.hpp
//  
//
//  Created by ShanshanLi on 7/24/15.
//
//

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

class dependency_graph
{
private:
    int    Nchl;   // number of channels
    int    Nb;      // predefined maximum of neighbors of one channel
    double diffusion_thres; // threshold of radius for interactions between channels
public:
    dependency_graph(int,int,double);
    boost::numeric::ublas::matrix<double>  channel_neighbors();
    boost::numeric::ublas::matrix<int>  Depend_graph();
};