//
//  Test_dependency_graph.cpp
//  
//
//  Created by ShanshanLi on 7/21/15.
//
//

#include "channel_coordinate.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/multi_array.hpp>
#include <boost/spirit/include/karma.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>


typedef boost::numeric::ublas::matrix<double> matrix;

//namespace karma = boost::spirit::karma;


/*template <typename T> boost::const_multi_array_ref<T, 2> make_view(boost::numeric::ublas::matrix<T> const& m) {
    return  boost::const_multi_array_ref<T,2> (
                                               &*m.data().begin(),
                                               boost::extents[m.size1()][m.size2()]
                                               );
}*/

void printPt(std::ostream &out, matrix p)			// print point
{
    
    for (int i = 0; i <p.size1(); i++) {
        for (int j=0; j<p.size2(); j++) {
             out<<std::setprecision(7)<<p(i,j)<<"\t";
        }
    out << "\n";
    }
}

int main(int argc, char **argv)
{
    double     radius=1;
    int        Nchl=2;
    channel_coordinate graph(Nchl,radius);
    matrix m= graph.coordinates();
    std::ofstream chal_dist;
    chal_dist.open("channel_coordinates.txt");
    //std::fill(m.data().begin(), m.data().end(), 1);
    
    
    //using namespace karma;
    //chal_dist<<format(auto_ % '\t' % eol,make_view(m)) << "\n";
    printPt(chal_dist, m);
    chal_dist.close();
    
    return 0;
}
