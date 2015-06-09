#include "diff_schemes.hpp"

// definition without declaration 'cause I'm too lazy
double eufunc(double t, double u){
    return u * cos(t);
}

double solu(double t){
    return exp(sin(t));
}

// a little bit stupid but well...
double pack_cos(double t, double u){
    return cos(t);
}

double warm(double t, double paul, double john){
    return john;
}

double gun(double t, double paul, double john){
    return -paul;
}


int main(){
    std::ofstream happiness; 
    imn::explicit_euler(eufunc, 1, M_PI/400., 0, 4.*M_PI, happiness, solu);
    imn::implicit_euler(eufunc, pack_cos, 1, M_PI/400., 0, 4.*M_PI, happiness, solu);
    imn::newton_euler(eufunc, pack_cos, 1, M_PI/400., 0, 4.*M_PI, happiness, 1e-6, solu);
    imn::richardson_euler(eufunc, 1, M_PI/400., 0, 4.*M_PI, happiness, solu);
    imn::autostep_euler(eufunc, 1, 5.*M_PI, 0, 4.*M_PI, happiness, 1e-5, 0.75, solu);
    imn::rk4(warm, gun, 0, 1, 0.1, 0, 4.*M_PI, happiness); 
}
