#define __LAB4__

#ifdef __LAB1__
#include "diff_schemes/diff_schemes.hpp"

// definition without declaration 'cause test purposes etc.
// differential equasion and its solution (for comparsion)
double eufunc(double t, double u){ return u * cos(t); }

double solu(double t){ return exp(sin(t)); }

// if your derivative don't need 'u', you have to do something like that
double pack_cos(double t, double u){ return cos(t); }

// rk4 functions with The Beatles reference
double warm(double t, double paul, double john){ return john; }

double gun(double t, double paul, double john){ return -paul; }
#endif

#ifdef __LAB4__
#include "vi_flow.hpp"
// vorticity function
double vor(double x, double y){ return -0.5 * (2 * y); }

// flux function
double flu(double x, double y){ return -0.5 * ((y*y*y)/3. - y*0.5*0.5); }
#endif

int main(){
    std::ofstream happiness;

    #ifdef __LAB1__
    imn::explicit_euler(eufunc, 1, M_PI/400., 0, 4.*M_PI, happiness, solu);
    imn::implicit_euler(eufunc, pack_cos, 1, M_PI/400., 0, 4.*M_PI, happiness, solu);
    imn::newton_euler(eufunc, pack_cos, 1, M_PI/400., 0, 4.*M_PI, happiness, 1e-6, solu);
    imn::richardson_euler(eufunc, 1, M_PI/400., 0, 4.*M_PI, happiness, solu);
    imn::autostep_euler(eufunc, 1, 5.*M_PI, 0, 4.*M_PI, happiness, 1e-5, 0.75, solu);
    imn::rk4(warm, gun, 0, 1, 0.1, 0, 4.*M_PI, happiness);
    #endif

    #ifdef __LAB4__
    imn::Point a,b,c,d;
    a = {-0.1, -0.5};
    b = {0.1, -0.5};
    c = {-0.1, 0};
    d = {0.1, 0};
    imn::Wall w1, w2, w3,w4;
    w1 = {a, c, imn::Wall::LEFT};
    w2 = {b, d, imn::Wall::RIGHT};
    w3 = {c, d, imn::Wall::UP};
    w4 = {a, b, imn::Wall::DOWN};
    std::vector<imn::Wall*> walls = {&w1, &w2, &w3, &w4};
    imn::poiseuille(-1, 2, -0.5, 0.5, 0.01, 0.01, flu, vor, happiness);
    imn::ns_obstacle(-1, 2, -0.5, 0.5, 0.01, 0.01, flu, vor, walls, happiness);
    #endif

    return 0;
}
