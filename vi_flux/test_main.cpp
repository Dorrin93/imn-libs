#include "vi_flux.hpp"

double vor(double x, double y){
    return -0.5 * (2 * y);
}

double flu(double x, double y){
    return -0.5 * ((y*y*y)/3. - y*0.5*0.5);
}

int main(){
    std::ofstream a;
    imn::poiseuille(-1, 2, -0.5, 0.5, 0.01, 0.01, flu, vor, a);
}
