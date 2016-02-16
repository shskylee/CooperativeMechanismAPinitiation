//
//  Created by ShanshanLi on 7/24/15.



#include "HeapTree.h"
#include "dependency_graph.hpp"
#include  "cdf.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <random>
#include <string>

using namespace boost::numeric::ublas;

/*
int    T_max=10000;          //millisecond
int    Steps_max=1000000;

// parameters for distribution of channels
int    Nchls=12008;          // number of channels  158 clusers * 2 rings per cluster * 38 channels per ring
double diffusion_thres=2.5;  // maximum distance for diffusion
int    Nbs=1000;             // maximum of neighbors of one channel
 */

//test algorithm{
const double T_max=1.0;             //millisecond
const int    Steps_max=1;
const int    Nchls=152;            // number of channels  2 clusers * 2 rings per cluster * 38 channels per ring
const double diffusion_thres=2.5;  // maximum distance for diffusion
const int    Nbs=100;              // maximum of neighbors of one channel
// }parameters for test


// parameter and functions for gating model of single sodium channel
const double reverse_voltage=-50.0;        // reverse sodium voltage millivolt
const int    Ns=3;                         // number of states per channel c(0),o(1),i(2)
//chosen reaction diagram in order: c->o (0), o->c (1), i->c (2), c->i(3), o->i (4)
const int    Nres=5;                      // number of reactions per channel
// single sodium channel gating parameters
const double V_A=-35.0;
const double k_A=6.0;
const double tau_A=0.1;
const double tau_I=0.5;
const double V_CI=80.0;
const double k_CI=4.0;
const double tau_CI=30.0;
const double pi = 3.1415926;


// ion diffusion parameters
const double Current_D=1.0;    // current released per channel
const double kappa=1.111e+3;        // inverse of the debye length (the debye length 9.0e-4)
const double epsilon=1.0;                // absolute permittivity of solution
const double Diff=1.0;                // diffusion constant in the solution


const double dt=1.0e-6;
const double err_bound=1.0e-6;

// print out function
void printPt(std::ostream &out, matrix<double> p);
// samples from exponential distribution random
double exp_dist(const double lambda);
// function of execuation of chosen reaction with least putative time
void exe_reaction(int x,int y,int z,int mu0, int n0);
// function to update open neighbors
void neighbor_update(int mu0, int n0, double t0, matrix<int> &matrix1, matrix<double> &matrix2, matrix<double> &matrix3);

int main()
{
    int    i,j,reaction;
    int    *column_sum;
    column_sum = new int[Nchls];
    int    counter=0;
    double tau_mu;
    double tau_prime;
    int    mu;
    int    channel;
    int    reactant;
    int    product;
    cdf    cumul_func;

    
    // first step : initialize;
    // (a) set initial states at t=0, dependency graph of reactions
    matrix<int>         system_state(Nchls*Ns,Steps_max);
    vector<double>      T(Steps_max);
    dependency_graph    graph(Nchls,Nbs,diffusion_thres);
    matrix<double>      Gmat=graph.Depend_graph();       // dependency graph
    matrix<double>      neighbors=graph.channel_neighbors();  // find out neighbors
    matrix<double>      neighbors_dist(Nchls,Nchls);
    matrix<double>      local_voltage(Nchls,Steps_max);
    matrix<int>         open_neighbors(Nchls,Nchls);   // neighbors that are open
    matrix<int>         open_neighbors_old(Nchls,Nchls);
    matrix<double>      open_times(Nchls,Nchls);       // times when neighbors become open
    matrix<double>      open_times_old(Nchls,Nchls);

    
    
    
    for (i=0; i<Nchls; i++) {
        system_state(i*Ns,0)=1;
        local_voltage(i,0)=reverse_voltage;
        for (j=0; j<Nbs; j++) {
            if (neighbors(i*Nbs+j,1)<=diffusion_thres)
                neighbors_dist(i,neighbors(i*Nbs+j,0))=neighbors(i*Nbs+j,1);
        }
    }
    


    //(b) calculate the propensity function a_i for all reactions
    int  Nreactions=Nchls*Nres;
    vector<double> propensity(Nreactions);
    
    for (i=0; i<Nchls; i++) {
        propensity(i*Nres+0)=system_state(i*Ns,0)*cumul_func.alpha_A(reverse_voltage);
        propensity(i*Nres+1)=system_state(i*Ns+1,0)*cumul_func.beta_A(reverse_voltage);
        propensity(i*Nres+2)=system_state(i*Ns+2,0)*cumul_func.alpha_CI(reverse_voltage);
        propensity(i*Nres+3)=system_state(i*Ns,0)*cumul_func.beta_CI(reverse_voltage);
        propensity(i*Nres+4)=system_state(i*Ns+1,0)/tau_I;
    }
    
    // (c) for each reaction i, generate a putative time, according to exponential distribution with parameter propensity a_i
    double * putative_times;
    putative_times= new double[Nreactions];
    
    for (i=0;i<Nreactions; i++) {
       putative_times[i]=exp_dist(propensity(i));
    }
    // (d) store the tau_i values in an indexed priority queue
    HeapTree<double> queue(putative_times,Nreactions,Nreactions);
    T(counter)=0;
    
    
    while (T(counter)<T_max&&counter<Steps_max) {

// second step: let mu be the reaction whose putative time tau_mu is least
       
        queue.Root_Node(tau_mu,mu);
        
// third & fourth steps: let tau be tau_mu, change the number of molecules to reflect execution of reaction mu, set t=tau_mu,update open neighbors, local voltage
       
        exe_reaction(reactant,product,channel,mu,Nres);
        system_state(channel*Ns+reactant,counter)-=1;
        system_state(channel*Ns+product,counter)+=1;
        counter++;
        T(counter)=tau_mu;
        neighbor_update(mu,Nres,tau_mu,open_neighbors,open_times,neighbors);
        
        for (i=0; i<Nchls; i++) {
            column_sum[i]=0;
            for (j=0; j<Nchls; j++) {
                column_sum[i] += open_neighbors_old(i,j);
            }
            if (!column_sum[i])
                local_voltage(i,counter)=cumul_func.local_voltage_update(i,T(counter-1),T(counter),local_voltage(i,counter-1),open_neighbors_old,open_times_old,neighbors_dist);
            else
                local_voltage(i,counter)=local_voltage(i,counter-1);
        }
        
        
// fifth step: (a)for each edge (mu,alpha) in the dependency graph Gmat, finally update the queue
        
        for (i=0; i<Nchls; i++) {
            for (j=0; j<Nres; j++) {
                reaction=i*Nres+j;
                if (!Gmat(mu,reaction)||!column_sum[i])
                {
                    cumul_func.set_values(i,j,system_state(reaction,counter-1),system_state(reaction,counter),T(counter-1),T(counter),putative_times[reaction],local_voltage(i,counter-1),local_voltage(i,counter),open_neighbors_old,open_neighbors,open_times_old,open_times,neighbors_dist);
                    tau_prime=cumul_func.inverse_cdf();
                    queue.Update(reaction,tau_prime);
                    
                }
            }
        }
     // neighbors and their open times
        neighbor_update(mu,Nres,tau_mu,open_neighbors_old,open_times_old,neighbors);
        
    // sixth step: go to step 2
}// end while
    
    
    
    delete[] putative_times;
    delete[] column_sum;
    
    return 0;
}

void printPt(std::ostream &out, matrix<double> p)
{
    
    for (int i = 0; i <p.size1(); i++) {
        for (int j=0; j<p.size2(); j++) {
            out<<std::setprecision(7)<<p(i,j)<<"\t";
        }
        out << "\n";
    }
}

double exp_dist(const double lambda)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::exponential_distribution<> distr(lambda);
    
    double rnd;
    
    return rnd=distr(gen);
    
}

void exe_reaction(int x,int y,int z,int mu0, int n0)
{
    int m=mu0%n0;
    z=(mu0-m)/n0;
    //chosen reaction diagram in order: c(0)->o(1) (0), o(1)->c(0) (1), i(2)->c(0) (2), c(0)->i(2) (3), o(1)->i(2) (4)
    switch (m) {
        case 0:
            x=0;
            y=1;
            break;
        case 1:
            x=1;
            y=0;
            break;
        case 2:
            x=2;
            y=0;
            break;
        case 3:
            x=0;
            y=2;
            break;
        case 4:
            x=1;
            y=2;
            break;
        default:
            break;
    }
}

void neighbor_update(int mu0, int n0, double t0, matrix<int> &matrix1, matrix<double> &matrix2, matrix<double> &matrix3)
{
    int m=mu0%n0;
    int z=(mu0-m)/n0;
    
    switch (m) {
        case 0:
            for (int j=0; j<Nbs; j++) {
                if (matrix3(z*Nbs+j,1)<=diffusion_thres) {
                    matrix1(matrix3(z*Nbs+j,0),z)=1;
                    matrix2(matrix3(z*Nbs+j,0),z)=t0;
                }
            }
            break;
        case 1:
           for (int j=0; j<Nbs; j++) {
                if (matrix3(z*Nbs+j,1)<=diffusion_thres) {
                    matrix1(matrix3(z*Nbs+j,0),z)=0;
                    matrix2(matrix3(z*Nbs+j,0),z)=0;
                }
            }
            break;
        case 4:
           for (int j=0; j<Nbs; j++) {
                if (matrix3(z*Nbs+j,1)<=diffusion_thres) {
                    matrix1(matrix3(z*Nbs+j,0),z)=0;
                    matrix2(matrix3(z*Nbs+j,0),z)=0;
                }
            }
            break;
        default:
            break;
    }
}






