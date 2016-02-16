//
//  cdf.h
//
//  cumulative distribution
//
//  Created by ShanshanLi on 8/10/15.
//
//


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <iostream>
#include <math.h>
#include <random>

extern const double V_A;
extern const double k_A;
extern const double tau_A;
extern const double tau_I;
extern const double V_CI;
extern const double k_CI;
extern const double tau_CI;
extern const double pi;
extern const double Current_D;
extern const double kappa;
extern const double epsilon;
extern const double Diff;
extern const double dt;
extern const double err_bound;



using namespace boost::numeric::ublas;

class cdf
{
public:
    // single sodium channel reaction propensity functions
    double alpha_A(double voltage);
    double beta_A(double voltage);
    double alpha_CI(double voltage);
    double beta_CI(double voltage);
    // local voltage increase functions due to opening channel itself (f1) and open neighbors(f2)
    double voltage_increase_f1(double t_end, double t_start);
    double voltage_increase_f2(double t_end, double t_start,double distance);
    // local voltage update
    double local_voltage_update(int a, double t_start, double t_end, double voltage_old, matrix<int> &m1, matrix<double> &m2, matrix<double> &m3);
    // set values
    void set_values(int a0, int a, double b, double c,double t0, double tm, double t1, double v0, double v1, matrix<int> &m1, matrix<int> &m2, matrix<double> &m3, matrix<double> &m4, matrix<double> &m5);
    // integrate
    double integrate(int state, double t_start, double t_end, double voltage_old,matrix<int> &m1, matrix<double> &m2, matrix<double> &m3);
    // probability distribution
    double prob_dist(int x, double tx);
    //cumulative distribution
    double cumul_dist(int x,double ty);
    // inverse of cumulative distribution
    double inverse_cdf(void);
    // destructor
    //~cdf(void);

private:
    
    int  state_old;
    int  state_new;
    double  t_old;                                   // time when (n-1)th event happened
    double t_middle;                                // time when  nth event happened
    double t_new;                             // time when putative time for considered reaction
    int    channel0;
    int    reaction0;                               // reaction we are considering
    double voltage_old0;                            // local voltage at time t_old around the channel0
    double voltage_new0;                            // local voltage at time t_new around the channel0
    matrix<int>   open_neighbors_old0;         // open neighbors of channel0 at t_old
    matrix<int>   open_neighbors_new0;         // open neighbors of channel0 at t_new
    matrix<double>  neighbors_distances0;     // neighbor distances of channel0
    matrix<double>  open_times_old0;          // open times of neighbors of channel0 at t_old
    matrix<double>  open_times_new0;          // open times of neighbors of channel0 at t_new
};


void cdf::set_values(int a0, int a, double b, double c,double t0, double tm, double t1, double v0, double v1, matrix<int> &m1, matrix<int> &m2, matrix<double> &m3, matrix<double> &m4, matrix<double> &m5)
{
    channel0=a0;
    reaction0=a;
    state_old=b;
    state_new=c;
    t_old=t0;
    t_middle=tm;
    t_new=t1;
    voltage_old0=v0;
    voltage_new0=v1;
    open_neighbors_old0=m1;
    open_neighbors_new0=m2;
    open_times_old0=m3;
    open_times_new0=m4;
    neighbors_distances0=m5;
}


// inverse of cumulative distribution
double cdf::inverse_cdf(void)
{
    // bisection method to find root
    // generate random lower bound and lower bound
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distr(0.0,1.0);

    double U,upper_bound,lower_bound,middle;
    double tau_new=0.0;
    double error;
    upper_bound=distr(gen);
    lower_bound=distr(gen);

    if (t_new>t_middle) {
        // generate new putative time for reaction rather a_mu
        U=(cumul_dist(0,t_new)-cumul_dist(0,t_middle))/(1.0-cumul_dist(0,t_middle));
        while ((cumul_dist(0,upper_bound)-U)<0.0) {
            upper_bound=distr(gen);
        }
        while ((cumul_dist(0,lower_bound)-U)>0.0) {
            lower_bound=distr(gen);
        }
        for (; ;) {
            middle=(upper_bound+lower_bound)/2.0;
            error=cumul_dist(0,middle)-U;
            if (error<err_bound) {
                tau_new=middle;
                return tau_new;
            }
            
            if ((cumul_dist(0,middle)-U)>0.0)
                upper_bound=middle;
            else
                lower_bound=middle;
            
        }
    }
    else
    {
     // generate new putative time for reaction a_mu
        U=distr(gen);
        while ((cumul_dist(1,upper_bound)-U)<0.0) {
            upper_bound=distr(gen);
        }
        while ((cumul_dist(1,lower_bound)-U)>0.0) {
            lower_bound=distr(gen);
        }
        for (; ;) {
            middle=(upper_bound+lower_bound)/2.0;
            error=cumul_dist(1,middle)-U;
            if (error<err_bound) {
                tau_new=middle;
                return tau_new;
            }
            
            if ((cumul_dist(1,middle)-U)>0.0)
                upper_bound=middle;
            else
                lower_bound=middle;
            
        }
    }

}


// cumulative probability distribution
double cdf::cumul_dist(int x, double ty)
{
    double integral0=0;
    int     steps0;
    if (x==0) {
        steps0=ceil((ty-t_old)/dt);
        for (int i=0; i<=steps0; i++) {
            integral0 += prob_dist(0,t_old+i*dt);
        }
    }
    else
    {
        steps0=ceil((ty-t_new)/dt);
        for (int i=0; i<=steps0; i++) {
            integral0 += prob_dist(1,t_old+i*dt);
        }
    }
    
    return integral0;
}

// probability distribution
double cdf::prob_dist(int x, double tx)
{
    double prob_d=0;
    
    //x=0 compute the cdf from t_n-1
    //x=1 compute the cdf from t_n
    
    if (x==0) {
        if (tx>t_middle) {
            prob_d=integrate(state_old,t_old,t_middle,voltage_old0,open_neighbors_old0,neighbors_distances0,open_times_old0);
            prob_d+=integrate(state_new,t_middle,tx,voltage_new0,open_neighbors_new0,neighbors_distances0,open_times_new0);
            prob_d = exp(-prob_d);
            prob_d *= state_new*local_voltage_update(channel0,t_middle,tx,voltage_new0,open_neighbors_new0,neighbors_distances0,open_times_new0);
        }
        else{
            prob_d=integrate(state_old,t_old,tx,voltage_old0,open_neighbors_old0,neighbors_distances0,open_times_old0);
            prob_d = exp(-prob_d);
            prob_d *= state_old*local_voltage_update(channel0,t_old,tx,voltage_old0,open_neighbors_old0,neighbors_distances0,open_times_old0);
        }
    }
    else
    {
        prob_d=integrate(state_new,t_middle,tx,voltage_new0,open_neighbors_new0,neighbors_distances0,open_times_new0);
        prob_d = exp(-prob_d);
        prob_d *= state_new*local_voltage_update(channel0,t_middle,tx,voltage_new0,open_neighbors_new0,neighbors_distances0,open_times_new0);
    }
    
       return prob_d;
}
        

// integration
double cdf::integrate(int state, double t_start, double t_end, double voltage_old,matrix<int> &m1, matrix<double> &m2, matrix<double> &m3)
{
    // state is either state_old or state_new
    double integral=0;
    int steps=ceil((t_end-t_start)/dt);
    double voltage_update=0;
    
    switch (reaction0) {
        case 0:
            for (int i=0; i<=steps; i++){
                voltage_update=local_voltage_update(channel0,t_start,t_start+i*dt,voltage_old,m1,m2,m3);
                integral += state*alpha_A(voltage_update)*dt;
            }
            break;
        case 1:
            for (int i=0; i<=steps; i++){
                voltage_update=local_voltage_update(channel0,t_start,t_start+i*dt,voltage_old,m1,m2,m3);
                integral += state*beta_A(voltage_update)*dt;
            }
            break;
        case 2:
            for (int i=0; i<=steps; i++){
                voltage_update=local_voltage_update(channel0,t_start,t_start+i*dt,voltage_old,m1,m2,m3);
                integral += state*alpha_CI(voltage_update)*dt;
            }
            break;
        case 3:
            for (int i=0; i<=steps; i++){
                voltage_update=local_voltage_update(channel0,t_start,t_start+i*dt,voltage_old,m1,m2,m3);
                integral += state*beta_CI(voltage_update)*dt;
            }
            break;
        case 4:
            for (int i=0; i<=steps; i++)
                 integral += state/tau_I*dt;
            
            break;
        default:
            break;
    }

    return integral;
    
}


double cdf::alpha_A(double voltage)
{
    double prop;
    prop=1/(1+exp(-(voltage-V_A)/k_A))/tau_A;
    return prop;
}


double cdf::beta_A(double voltage)
{
    double prop;
    prop=1/(1+exp((voltage-V_A)/k_A))/tau_A;
    return prop;
}

double cdf::alpha_CI(double voltage)
{
    double prop;
    prop=(1+exp(-(voltage-V_CI)/k_CI))/tau_CI;
    return prop;
}


double cdf::beta_CI(double voltage)
{
    double prop;
    prop=(1+exp((voltage-V_CI)/k_CI))/tau_CI;
    return prop;
}
        

double cdf::voltage_increase_f1(double t_end, double t_start)
{
    // voltage_increase_f1: function for local voltage increase due to the opening of channel itself
    double V_increase=0;
    
    V_increase=Current_D/(4*Diff*pi*epsilon*kappa)*(1-exp(Diff*(t_end-t_start)*pow(kappa,2.0))*erfc(sqrt(Diff*(t_end-t_start))*kappa));
    
    return V_increase;
}


double cdf::voltage_increase_f2(double t_end, double t_start,double distance)
{
    // voltage_increase_f2: function for local voltage increase due to the opening of neighbors
    double V_increase=0;
    
    V_increase=(exp(1.0)-2)*Current_D/(4*Diff*exp(1.0)*pi*epsilon*pow(kappa,2.0)*distance)*erfc(distance/(2*sqrt(Diff*(t_end-t_start))));
    
    return V_increase;
}
        

double cdf::local_voltage_update(int a, double t_start, double t_end, double voltage_old,matrix<int> &m1, matrix<double> &m2, matrix<double> &m3)

{
    
    // this subroutine we calculate increase of local voltage a the ath channel per step
   
    // this increase is affected by how many neighbors are open, the time difference from their open times, distances
   
    // column1 is the open neighbors of ath channel at t_start
    // column2 is the open times of neighbor channels at to t_start
    // column3 is the neighbors distances of ath channel
    
    // voltage_increase_f1: function for local voltage increase due to the opening of channel itself
    // voltage_increase_f2: function for local voltage increase due to the opening of neighbors
    
    double local_voltage=voltage_old;
    
    matrix_row<matrix<int> > vec1(m1,a);
    matrix_row<matrix<double> > vec2(m2,a);
    matrix_row<matrix<double> > vec3(m3,a);
    
    for (int i=0; i<vec1.size(); i++) {
        if (i==a) {
            local_voltage += voltage_increase_f1(t_end,vec2(i))-voltage_increase_f1(t_start,vec2(i)); // contributions from opening of itself
        }
        else if (!vec1(i)){
            local_voltage += voltage_increase_f2(t_end,vec2(i),vec3(i))-voltage_increase_f2(t_start,vec2(i),vec3(i)); // contributions from opening neighbors
        }
        else
            local_voltage += 0; // no contributions
    }
    
    return local_voltage;
}






