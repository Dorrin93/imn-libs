#include "vi_flux.hpp"


double imn::vortHelper(double **f, double **v, int x, int y)
{
    // 0.0625 is 1/16
    return 0.25 * (v[x+1][y] + v[x-1][y] + v[x][y-1] + v[x][y+1]) -
        0.0625 * ( (f[x][y+1] - f[x][y-1]) * (v[x+1][y] - v[x-1][y]) -
        (f[x+1][y] - f[x-1][y]) * (v[x][y+1] - v[x][y-1]) );
}


double imn::fluxHelper(double **f, double **v, int x, int y,
        double dx, double dy)
{
    return (f[x+1][y] + f[x-1][y] + f[x][y+1] + f[x][y-1] - v[x][y] * 
            dx * dy) * 0.25;
}


bool imn::inside(double x, double y, const std::vector<Wall*> &walls)
{
    for(std::vector<Wall*>::const_iterator it = walls.begin();
            it != walls.end(); ++it)
    {
        switch((*it)->type){
            case UP: 
                if(x >= (*it)->A.x && x <= (*it)->B.x && y <= (*it)->A.y) 
                    return true; 
                break;
            case DOWN: 
                if(x >= (*it)->A.x && x <= (*it)->B.x && y >= (*it)->A.y) 
                    return true;
                break;
            case LEFT:
                if(y >= (*it)->A.y && y <= (*it)->B.y && x >= (*it)->A.x)
                    return true;
                break;
            case RIGHT:
                if(y >= (*it)->A.y && y <= (*it)->B.y && x <= (*it)->A.x)
                    return true;
                break;
        }
    }
    return false;
}


/////////////////////////////////////////////////////////////////////////////
void imn::poiseuille(double xmin, double xmax, double ymin,
        double ymax, double dx, double dy, imn::func flux,
        imn::func vort, std::ofstream &file, double prec, int min_iter)
{
    // 'gettin grid size
    int x_size = static_cast<int>(round( (xmax - xmin) / dx ));
    int y_size = static_cast<int>(round( (ymax - ymin) / dy ));

    // gird allocation
    double** f_grid = new double*[x_size+1];
    double** v_grid = new double*[x_size+1];
    for(int i = 0; i <= x_size; ++i){
        f_grid[i] = new double[y_size+1];
        v_grid[i] = new double[y_size+1];
    }

    // initial conditions
    // inside
    for(int i = 0; i <= x_size; ++i){
        for(int j = 0; j <= y_size; ++j){
            f_grid[i][j] = 0;
            v_grid[i][j] = 0;
        }
    }
    // up and down
    for(int i = 0; i <= x_size; ++i){
        f_grid[i][0] = flux(xmin + i * dx, ymin);
        f_grid[i][y_size] = flux(xmin + i * dx, ymin + y_size * dy);
        v_grid[i][0] = vort(xmin + i * dx, ymin);
        v_grid[i][y_size] = vort(xmin + i * dx, ymin + y_size * dy);
    }
    // left and right
    for(int j = 0; j <= y_size; ++j){
        f_grid[0][j] = flux(xmin, ymin + j * dy);
        f_grid[x_size][j] = flux(xmin + x_size * dx, ymin + j * dy);
        v_grid[0][j] = vort(xmin, ymin + j * dy);
        v_grid[x_size][j] = vort(xmin + x_size * dx, ymin + j * dy);
    }

    // some special values - we check condition in 2/3 x and 0.5 y point
    int p_x = static_cast<int>(round( (2./3.) * x_size ));
    int p_y = y_size >> 1; // so fancy...
    double prev_f = 0;
    double prev_v = 0;
    int iter = 0;

    while(iter < min_iter || fabs(f_grid[p_x][p_y] - prev_f) > prec ||
            fabs(v_grid[p_x][p_y] - prev_v) > prec){
        
        ++iter;
        prev_f = f_grid[p_x][p_y];
        prev_v = v_grid[p_x][p_y];

        // I've heard you like loops, so I've put the loop into the loop
        // into the loop, so you can iterate while you iterate while
        // you acually calculate, so... yeah that's n^3 nightmare.
        
        // vorticity comes first
        for(int i = 1; i < x_size; ++i){
            for(int j = 1; j < y_size; ++j){
                v_grid[i][j] = imn::vortHelper(f_grid, v_grid, i, j);
            }
        }
        // flux is next (they cannot be done in same loop)
        for(int i = 1; i < x_size; ++i){
            for(int j = 1; j < y_size; ++j){
                f_grid[i][j] = imn::fluxHelper(f_grid, v_grid, i, j , dx, dy);
            }
        }
    }

    // checking if file is was open and eventually opening it
    bool was_open = file.is_open();
    if(!was_open)
        file.open("poiseuille.dat");

    file << "# x y u(x,y) v(x,y)" << std::endl;

    for(int i = 0; i < x_size; ++i){
        for(int j = 0; j < y_size; ++j){
            file << xmin + i*dx << " " << ymin + j*dy << " "
                 << (f_grid[i][j+1] - f_grid[i][j]) / dy  << " " 
                 << -(f_grid[i+1][j] - f_grid[i][j]) / dx 
                 << std::endl;
        }
        file << std::endl;
    }

    //if file wasn't open closing it
    if(!was_open)
        file.close();


    // grid deallocation
    for(int i = 0; i <= x_size; ++i){
        delete[] f_grid[i];
        delete[] v_grid[i];
    }
    delete[] f_grid;
    delete[] v_grid;

}

