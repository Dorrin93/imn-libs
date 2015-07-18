/**
 * @file   grid.hpp
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   July, 2015
 * @brief  Implemenation of 2D obstackle grid class.
 */
#include "ogrid.hpp"


void imn::OGrid::apply_to_obstackle(func2d function) noexcept
{
    if(obstacle){

        auto i = 0u;
        auto j = 0u;

        for(auto& it : mtx_){

            if(obstacle(xpos(i), ypos(j)))
                it = function(xpos(i), ypos(j));

            increm(i, j);
        }
    }
    else{
        std::cerr << "Obstackle funtion not given" << std::endl;
    }
}

void imn::OGrid::apply_not_to_obstackle(func2d function) noexcept
{
    if(obstacle){

        for(auto i = 1u; i < x_size()-1; ++i)
            for(auto j = 1u; j < y_size()-1; ++j) {

                if(!obstacle(xpos(i), ypos(j)))
                    (*this)(i, j) = function(xpos(i), ypos(j));
            }
    }
    else{
        std::cerr << "Obstackle funtion not given" << std::endl;
    }
}

void imn::OGrid::auto_von_Neumann(VonNeumann type)
{
    auto func = static_cast<vNfunc>(nullptr);

    switch(type){
        case VonNeumann::pFlow:
            func = pFlowFunc;
            break;

        case VonNeumann::viFlow:
            func = viFlowFunc;
            break;
    }

    if(obstacle){

        for(auto i = 1u; i < x_size()-1; ++i){
            for(auto j = 1u; j < y_size()-1; ++j){

                if( obstacle(xpos(i), ypos(j)) )
                    (*this)(i, j) = func(i, j, type(i, j));

            }
        }

    }
    else{
        std::cerr << "Obstackle funtion not given" << std::endl;
    }

}


void imn::OGrid::von_Neumann_cond(imn::func2d function, imn::OGrid::Ptype pointType)
{
    if(obstacle){

        for(auto i = 1u; i < x_size()-1; ++i){
            for(auto j = 1u; j < y_size()-1; ++j){

                if( obstacle(xpos(i), ypos(j)) && type(i, j) == pointType)
                    (*this)(i, j) = function(xpos(i), ypos(j));

            }
        }

    }
    else{
        std::cerr << "Obstackle funtion not given" << std::endl;
    }
}

imn::OGrid::Ptype imn::OGrid::type(unsigned i, unsigned j) const noexcept
{
    auto left = obstacle(xpos(i-1), ypos(j));
    auto right = obstacle(xpos(i+1), ypos(j));
    auto up = obstacle(xpos(i), ypos(j+1));
    auto down = obstacle(xpos(i), ypos(j-1));

    if(left && right && up && down) return Ptype::in;
    if(!left && !up)                return Ptype::lu;
    if(!left && !down)              return Ptype::ld;
    if(!right && !up)               return Ptype::ru;
    if(!right && !down)             return Ptype::rd;
    if(!left)                       return Ptype::left;
    if(!right)                      return Ptype::right;
    if(!up)                         return Ptype::up;
    if(!down)                       return Ptype::down;

    return Ptype::in;
}

double imn::OGrid::pFlowFunc(unsigned i, unsigned j, imn::OGrid::Ptype type) const noexcept
{
    switch (type){
        case Ptype::left : return (*this)(i-1, j);
        case Ptype::right: return (*this)(i+1, j);
        case Ptype::up   : return (*this)(i, j+1);
        case Ptype::down : return (*this)(i, j-1);
        case Ptype::lu   : return 0.5 * ( (*this)(i-1, j) + (*this)(i, j+1) );
        case Ptype::ld   : return 0.5 * ( (*this)(i-1, j) + (*this)(i, j-1) );
        case Ptype::ru   : return 0.5 * ( (*this)(i+1, j) + (*this)(i, j+1) );
        case Ptype::rd   : return 0.5 * ( (*this)(i+1, j) + (*this)(i, j-1) );
        case Ptype::in   : return (*this)(i, j);
    }
}

double imn::OGrid::viFlowFunc(unsigned i, unsigned j, imn::OGrid::Ptype type) const noexcept
{
    //TODO
    return 0;
}


