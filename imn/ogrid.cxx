/**
 * @file   grid.hpp
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   July, 2015
 * @brief  Implemenation of 2D obstackle grid class.
 */
#include "ogrid.hpp"

imn::OGrid::OGrid(double x_min, double x_max, double y_min, double y_max, double dx, double dy, bool widen) :
                Grid(x_min, x_max, y_min, y_max, dx, dy, nullptr, widen) {}

void imn::OGrid::set_obstackle(ofunc obs){

    obstacle = obs;
    auto p = 1e-8;
    obstacle_p_m = [obs, p](double x, double y){ return obs(x-p, y-p) ||
                                                             obs(x+p, y-p) ||
                                                             obs(x-p, y+p) ||
                                                             obs(x+p, y+p); };
    auto set = false;

    for(auto i = 0u; i < x_size(); ++i){

        if(!set) {

            if (obstacle(xpos(i), ypos(0)) || obstacle(xpos(i), ypos(1))) {
                obs_type_ = Obstype::DOWN;
                set = true;
            }

            if (obstacle(xpos(i), ypos(y_size() - 1)) || obstacle(xpos(i), ypos(y_size() - 2))) {
                obs_type_ = Obstype::UP;
                set = true;
            }
        }

        for(auto j = 0u; j < y_size(); ++j){
            if(obstacle_p_m(xpos(i), ypos(j)))
                obs_in_.emplace_back(i, j, point_type(i, j));
            else
                obs_out_.emplace_back(i, j);
        }

    }
    /*
    for(auto& it : obs_in_){
        std::cout << std::get<0>(it) << " " << std::get<1>(it) << " ";
        switch (std::get<2>(it)){
            case Ptype::LEFT : std::cout << "left"; break;
            case Ptype::RIGHT: std::cout << "right"; break;
            case Ptype::UP   : std::cout << "up"; break;
            case Ptype::DOWN : std::cout << "down"; break;
            case Ptype::LU   : std::cout << "lu"; break;
            case Ptype::LD   : std::cout << "ld"; break;
            case Ptype::RU   : std::cout << "ru"; break;
            case Ptype::RD   : std::cout << "rd"; break;
            case Ptype::IN   : std::cout << "in"; break;
        }
        std::cout << std::endl;
    }
     */

    if(!set)
        obs_type_ = Obstype::MID;
}

void imn::OGrid::apply_to_obstackle(func2d function) noexcept
{

    if(obstacle){

        for(const auto& it : obs_in_){
            (*this)(std::get<0>(it), std::get<1>(it)) = function(xpos(std::get<0>(it)), ypos(std::get<1>(it)));
        }

    }
    else{
        std::cerr << "Obstackle funtion not given" << std::endl;
    }
}

void imn::OGrid::apply_not_to_obstackle(func2d function) noexcept
{
    if(obstacle){

        for(const auto& it : obs_out_){
            if(it.first != 0 && it.first != x_size()-1 && it.second != 0 && it.second != y_size()-1)
                (*this)(it.first, it.second) = function(xpos(it.first), ypos(it.second));
        }
    }
    else{
        std::cerr << "Obstackle funtion not given" << std::endl;
    }
}

void imn::OGrid::apply_index_not_to_obstackle(gridFunc function) noexcept
{
    if(obstacle){

        for(const auto& it : obs_out_){
            if(it.first != 0 && it.first != x_size()-1 && it.second != 0 && it.second != y_size()-1)
                (*this)(it.first, it.second) = function(it.first, it.second);
        }
    }
    else{
        std::cerr << "Obstackle funtion not given" << std::endl;
    }
}

void imn::OGrid::auto_von_Neumann(Condtype type, const OGrid *other)
{
    auto func = static_cast<vNfunc>(nullptr);

    switch(type){
        case Condtype::PFLOW:
            func = [&,this](unsigned i, unsigned j, Ptype t, const OGrid* o){ return this->pFlowFunc(i, j, t, o); };
            break;

        case Condtype::VIFLOW:
            func = [&,this](unsigned i, unsigned j, Ptype t, const OGrid* o){ return this->viFlowFunc(i, j, t, o); };
            break;
    }

    if(obstacle){

        for(const auto& it : obs_in_){
            (*this)(std::get<0>(it), std::get<1>(it)) = func(std::get<0>(it), std::get<1>(it), std::get<2>(it), other);
        }

        switch(type){
            case Condtype::PFLOW:
                for(auto i = 0u; i < x_size(); ++i){
                    (*this)(i, 0) = (*this)(i, 1);
                    (*this)(i, y_size()-1) = (*this)(i, y_size()-2);
                }
                break;
            case Condtype::VIFLOW:
                for(auto i = 0u; i < x_size(); ++i){
                    (*this)(i, 0) = 2 * ((*other)(i, 1) - (*other)(i, 0)) / (dx() * dy());
                    (*this)(i, y_size()-1) = 2 * ((*other)(i, y_size()-2) - (*other)(i, y_size()-1)) / (dx() * dy());
                }
                break;
        }

    }
    else{
        std::cerr << "Obstackle funtion not given" << std::endl;
    }

}

void imn::OGrid::von_Neumann_cond(imn::func2d function, imn::OGrid::Ptype pointType)
{
    if(obstacle){

        for(const auto& it : obs_in_){
            if(std::get<2>(it) == pointType)
                (*this)(std::get<0>(it), std::get<1>(it)) = function( xpos(std::get<0>(it)), ypos(std::get<1>(it)) );
        }

    }
    else{
        std::cerr << "Obstackle funtion not given" << std::endl;
    }
}

imn::OGrid::Ptype imn::OGrid::point_type(unsigned i, unsigned j) const noexcept
{
    auto left = obstacle_p_m(xpos(i-1), ypos(j));
    auto right = obstacle_p_m(xpos(i+1), ypos(j));
    auto up = obstacle_p_m(xpos(i), ypos(j+1));
    auto down = obstacle_p_m(xpos(i), ypos(j-1));

    if(left && right && up && down) return Ptype::IN;
    if(!left && !up)                return Ptype::LU;
    if(!left && !down)              return Ptype::LD;
    if(!right && !up)               return Ptype::RU;
    if(!right && !down)             return Ptype::RD;
    if(!left)                       return Ptype::LEFT;
    if(!right)                      return Ptype::RIGHT;
    if(!up)                         return Ptype::UP;
    if(!down)                       return Ptype::DOWN;

    return Ptype::IN;
}

double imn::OGrid::pFlowFunc(unsigned i, unsigned j, Ptype type, const OGrid *other) const noexcept
{
    switch (type){
        case Ptype::LEFT : return (*this)(i-1, j);
        case Ptype::RIGHT: return (*this)(i+1, j);
        case Ptype::UP   : return (*this)(i, j+1);
        case Ptype::DOWN : return (*this)(i, j-1);
        case Ptype::LU   : return 0.5 * ( (*this)(i-1, j) + (*this)(i, j+1) );
        case Ptype::LD   : return 0.5 * ( (*this)(i-1, j) + (*this)(i, j-1) );
        case Ptype::RU   : return 0.5 * ( (*this)(i+1, j) + (*this)(i, j+1) );
        case Ptype::RD   : return 0.5 * ( (*this)(i+1, j) + (*this)(i, j-1) );
        case Ptype::IN   : return (*this)(i, j);
    }
}

#define __TEST__

double imn::OGrid::viFlowFunc(unsigned i, unsigned j, Ptype type, const OGrid *other) const noexcept
{
    switch (type){
        case Ptype::IN   : return (*this)(i, j);
        case Ptype::LEFT : return 2 * ((*other)(i-1, j) - (*other)(i, j)) / (dx() * dy());
        case Ptype::RIGHT: return 2 * ((*other)(i+1, j) - (*other)(i, j)) / (dx() * dy());
        case Ptype::UP   : return 2 * ((*other)(i, j+1) - (*other)(i, j)) / (dx() * dy());
        case Ptype::DOWN : return 2 * ((*other)(i, j-1) - (*other)(i, j)) / (dx() * dy());

        #ifndef __TEST__
        case Ptype::LU   : return 0.5 * (viFlowFunc(i, j, Ptype::LEFT, other) + viFlowFunc(i, j, Ptype::UP, other));
        case Ptype::LD   : return 0.5 * (viFlowFunc(i, j, Ptype::LEFT, other) + viFlowFunc(i, j, Ptype::DOWN, other));
        case Ptype::RU   : return 0.5 * (viFlowFunc(i, j, Ptype::RIGHT, other) + viFlowFunc(i, j, Ptype::UP, other));
        case Ptype::RD   : return 0.5 * (viFlowFunc(i, j, Ptype::RIGHT, other) + viFlowFunc(i, j, Ptype::DOWN, other));
        #endif

        // simpler option for testing purposes
        #ifdef __TEST__
        case Ptype::LU   : return viFlowFunc(i, j, Ptype::UP, other);
        case Ptype::RU   : return viFlowFunc(i, j, Ptype::UP, other);
        case Ptype::LD   : return viFlowFunc(i, j, Ptype::DOWN, other);
        case Ptype::RD   : return viFlowFunc(i, j, Ptype::DOWN, other);
        #endif

    }
}
