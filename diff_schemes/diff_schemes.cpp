#include "diff_schemes.hpp"

using namespace imn;
/////////////////////////////////////////////////////////////////////////////////////
void explicit_euler(input inp_func, const double icondition, const double dt, 
        const double begin, const double end, std::ofstream &file, solution sol_func)
{
    // checking difference between solution(begin) and condition,
    // if more then about 0 user probably uses bad input
    if(sol_func && fabs(sol_func(begin) - icondition) > 1e-9){
        std::cerr << "Difference between solution function and initial \
            condition is " << sol_func(begin) - icondition << std::endl;
        std::cerr << "It should be about 0, so input could be incorrect."
            << std::endl;
    }

    // checking if file is was open and eventually opening it
    bool was_open = file.is_open();
    if(!was_open)
        file.open("explicit_euler.dat");

    // writing initial value
    file << begin << " " << icondition;
    if(sol_func) 
        file << " " << sol_func(begin) << " " << sol_func(begin) - icondition;
    file << std::endl;

    // initial variable value
    double u_n = icondition;

    for(double i = begin + dt; i <= end; i += dt){
        
        // value calculation in this moment
        u_n = u_n + dt * inp_func(i - dt, u_n);

        //writing this moment values
        file << i << " " << u_n;
        if(sol_func)
            file << " " << sol_func(i) << " " << sol_func(i) - u_n;
        file << std::endl;

    }

    //if file wasn't open closing it
    if(!was_open)
        file.close();

}

/////////////////////////////////////////////////////////////////////////////////////
void implicit_euler(input inp_func, input inp_deriv, const double icondition, 
        const double dt, const double begin, const double end, std::ofstream &file, 
        solution sol_func)
{
    // checking difference between solution(begin) and condition,
    // if more then about 0 user probably uses bad input
    if(sol_func && fabs(sol_func(begin) - icondition) > 1e-9){
        std::cerr << "Difference between solution function and initial \
            condition is " << sol_func(begin) - icondition << std::endl;
        std::cerr << "It should be about 0, so input could be incorrect."
            << std::endl;
    }

    //checking if file is was open and eventually opening it
    bool was_open = file.is_open();
    if(!was_open)
        file.open("implicit_euler.dat");

    // writing initial value
    file << begin << " " << icondition;
    if(sol_func) 
        file << " " << sol_func(begin) << " " << sol_func(begin) - icondition;
    file << std::endl;

    // initial variable value
    double u_n = icondition;

    for(double i = begin + dt; i <= end; i += dt){

        // value calculation
        u_n = u_n / (1. - dt * inp_deriv(i, u_n));

        // writing values to file
        file << i << " " << u_n;
        if(sol_func)
            file << " " << sol_func(i) << " " << sol_func(i) - u_n;
        file << std::endl;

    }

    //if file wasn't open closing it
    if(!was_open)
        file.close();

}

/////////////////////////////////////////////////////////////////////////////////////
void newton_euler(input inp_func, input inp_deriv, const double icondition,
        const double dt, const double begin, const double end, 
        std::ofstream &file, const double stop, solution sol_func)
{
    // checking difference between solution(begin) and condition,
    // if more then about 0 user probably uses bad input
    if(sol_func && fabs(sol_func(begin) - icondition) > 1e-9){
        std::cerr << "Difference between solution function and initial \
            condition is " << sol_func(begin) - icondition << std::endl;
        std::cerr << "It should be about 0, so input could be incorrect."
            << std::endl;
    }

    //checking if file is was open and eventually opening it
    bool was_open = file.is_open();
    if(!was_open)
        file.open("newton_euler.dat");

    // writing initial value
    file << begin << " " << 0 << " " << icondition;
    if(sol_func) 
        file << " " << sol_func(begin) << " " << sol_func(begin) - icondition;
    file << std::endl;

    // initial variable value
    double u_n = icondition;
    // value from external iteration before
    double u_pn = u_n;
    // value from internal iteration before
    double u_wn = u_n;
    // iteration counter
    int iter = 0;

    for(double i = begin + dt; i <= end; i += dt, iter = 0, u_pn = u_n){

        //u_n = u_wn - (u_wn - u_pn - dt * inp_func(i, u_wn)) / (1. - dt * inp_deriv(i, u_wn));
        // iteration
        while(stop < fabs(u_n - u_wn) || !iter){
            u_wn = u_n;

            //value calculation
            u_n = u_wn - (u_wn - u_pn - dt * inp_func(i, u_wn)) / (1. - dt * inp_deriv(i, u_wn));

            ++iter; 
        }

        // writing to file
        file << i << " " << iter << " " << u_n;
        if(sol_func)
            file << " " << sol_func(i) << " " << sol_func(i) - u_n;
        file << std::endl;

    }

    //if file wasn't open closing it
    if(!was_open)
        file.close();

}

/////////////////////////////////////////////////////////////////////////////////////
void richardson_euler(input inp_func, const double icondition, const double dt,
        const double begin, const double end, std::ofstream &file, solution sol_func)
{
    // checking difference between solution(begin) and condition,
    // if more then about 0 user probably uses bad input
    if(sol_func && fabs(sol_func(begin) - icondition) > 1e-9){
        std::cerr << "Difference between solution function and initial \
            condition is " << sol_func(begin) - icondition << std::endl;
        std::cerr << "It should be about 0, so input could be incorrect."
            << std::endl;
    }

    // checking if file is was open and eventually opening it
    bool was_open = file.is_open();
    if(!was_open)
        file.open("richardson_euler.dat");

    // writing initial value
    file << begin << " " << icondition;
    if(sol_func) 
        file << " " << sol_func(begin) << " " << sol_func(begin) - icondition;
    file << std::endl;

    // initial variable value
    double u_n = icondition;

    for(double i = begin + 2*dt; i <= end; i += 2*dt){

        // first short step
        double u_1 = u_n + dt * inp_func(i - 2 * dt, u_n);
        // second short step
        double u_2 = u_1 + dt * inp_func(i - dt, u_1);
        // long step
        double u_l = u_n + 2 * dt * inp_func(i - 2 * dt, u_n);
        // value
        u_n = 2 * u_2 - u_l;

        //writing this moment values
        file << i << " " << u_n;
        if(sol_func)
            file << " " << sol_func(i) << " " << sol_func(i) - u_n;
        file << std::endl;

    }

    //if file wasn't open closing it
    if(!was_open)
        file.close();

}

/////////////////////////////////////////////////////////////////////////////////////
void autostep_euler(input inp_func, const double icondition, const double init_dt,
        const double begin, const double end, std::ofstream &file,
        const double tolerance, const double S, solution sol_func)
{
    // checking S parameter condition
    if(S >= 1.){
        std::cerr << "S parameter must be S < 1 !" << std::endl;
        return;
    }

    // checking difference between solution(begin) and condition,
    // if more then about 0 user probably uses bad input
    if(sol_func && fabs(sol_func(begin) - icondition) > 1e-9){
        std::cerr << "Difference between solution function and initial \
            condition is " << sol_func(begin) - icondition << std::endl;
        std::cerr << "It should be about 0, so input could be incorrect."
            << std::endl;
    }

    // checking if file is was open and eventually opening it
    bool was_open = file.is_open();
    if(!was_open)
        file.open("autostep_euler.dat");

    // writing initial value
    file << begin << " " << init_dt << " " << icondition;
    if(sol_func) 
        file << " " << sol_func(begin) << " " << sol_func(begin) - icondition;
    file << std::endl;

    // initial variable value
    double u_n = icondition;
    // initial step value
    double dt = init_dt;
    // safety counter
    int cnt = 0;

    for(double i = begin + 2*dt; i <= end || cnt < 5; i += 2*dt, ++cnt){

        // first short step
        double u_1 = u_n + dt * inp_func(i - 2 * dt, u_n);
        // second short step and perhabs value
        double u_2 = u_1 + dt * inp_func(i - dt, u_1);
        // long step
        double u_l = u_n + 2 * dt * inp_func(i - 2 * dt, u_n);

        if(fabs(u_2 - u_l) < tolerance){
            u_n = u_2;

            // writing this moment values
            file << i << " " << dt << " " << u_n;
            if(sol_func)
                file << " " << sol_func(i) << " " << sol_func(i) - u_n;
            file << std::endl;
        }
        else{
            // getting back in time :3
            i -= 2 * dt;
        }

        // calculating step
        dt = sqrt( (S * tolerance) / fabs(u_2 - u_l) ) * dt;

    }

    //if file wasn't open closing it
    if(!was_open)
        file.close();

}


/////////////////////////////////////////////////////////////////////////////////////
void rk4(rkfunc equasion1, rkfunc equasion2, const double cond_x, const double cond_y,
        const double dt, const double begin, const double end, std::ofstream &file)
{
    // checking if file is was open and eventually opening it
    bool was_open = file.is_open();
    if(!was_open)
        file.open("rk4.dat");

    // writing initial values
    file << begin << " " << cond_x << " " << cond_y << std::endl;

    // initial values
    double x = cond_x;
    double y = cond_y;
    

    for(double i = begin + dt; i <= end; i += dt){

        // some calculations tamdamdam
        double k1x = equasion1(i - dt, x, y);
        double k1y = equasion2(i - dt, x, y);

        double k2x = equasion1(i - dt + dt / 2, x + dt * k1x / 2, y + dt * k1y / 2);
        double k2y = equasion2(i - dt + dt / 2, x + dt * k1x / 2, y + dt * k1y / 2);

        double k3x = equasion1(i - dt + dt / 2, x + dt * k2x / 2, y + dt * k2y / 2);
        double k3y = equasion2(i - dt + dt / 2, x + dt * k2x / 2, y + dt * k2y / 2);

        double k4x = equasion1(i, x + dt * k3x, y + dt * k3y);
        double k4y = equasion2(i, x + dt * k3x, y + dt * k3y);

        // and here we are...
        x = x + (dt / 6.) * (k1x + 2 * k2x + 2 * k3x + k4x);
        y = y + (dt / 6.) * (k1y + 2 * k2y + 2 * k3y + k4y);

        file << i << " " << x << " " << y << std::endl;

    }

    //if file wasn't open closing it
    if(!was_open)
        file.close();
}

