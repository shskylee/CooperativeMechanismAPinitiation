//
//  Test_dependency_graph.cpp
//  
//
//  Created by ShanshanLi on 7/24/15.
//
//  number of clusters: 158
//  number of channels per cluster: 76
//  number of rings per cluster: 2
//  number of channels per ring: 38
//  angular distance between channels in one cluster :2*pi/38
//  distance between clusters: 190 nm 0.19 micrometers

// the diameter of axon initial segment: 1 micrometer
// the length of axon initial segment: 30 micrometers


#include "dependency_graph.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/multi_array.hpp>
#include <boost/spirit/include/karma.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

//namespace karma = boost::spirit::karma;

/*template <typename T>
boost::const_multi_array_ref<T, 2> make_view(boost::numeric::ublas::matrix<T> const& m) {
    return  boost::const_multi_array_ref<T,2> (
                                               &*m.data().begin(),
                                               boost::extents[m.size1()][m.size2()]
                                               );
}*/

typedef boost::numeric::ublas::matrix<double> matrix;

void printPt(std::ostream &out, matrix p)			// print point
{
    
    for (int i = 0; i <p.size1(); i++) {
        for (int j=0; j<p.size2(); j++) {
            out<<std::setprecision(7)<<p(i,j)<<"\t";
        }
        out << "\n";
    }
}



int main()
{
    int    Nchls=152;
    double radius_thres=2.5;
    int    Nbs=100;  // maximum of neighbors of one channel
    dependency_graph graph(Nchls,Nbs,radius_thres);
    boost::numeric::ublas::matrix<double> m= graph.Depend_graph();
    
   /* using namespace karma;
    std::cout << format(auto_ % '\t' % eol, make_view(m)) << "\n";*/
    printPt(std::cout,m);
    //std::cout<<m.size2()<<std::endl;
    return 0;
}