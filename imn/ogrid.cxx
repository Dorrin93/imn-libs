/**
 * @file   grid.hpp
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   July, 2015
 * @brief  Implemenation of 2D obstackle grid class.
 */
#include "ogrid.hpp"


void imn::OGrid::set_obstackle(ofunc obs){

    obstacle = obs;

    for(auto i = 0u; i < x_size(); ++i){

        if(obstacle(xpos(i), ypos(0)) || obstacle(xpos(i), ypos(1))){
            obs_type_ = Obstype::DOWN;
            return;
        }

        if(obstacle(xpos(i), ypos(y_size()-1)) || obstacle(xpos(i), ypos(y_size()-2))) {
            obs_type_ = Obstype::UP;
            return;
        }

    }

    obs_type_ = Obstype::MID;
}

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

void imn::OGrid::apply_index_not_to_obstackle(gridFunc function) noexcept
{
    if(obstacle){

        for(auto i = 1u; i < x_size()-1; ++i)
            for(auto j = 1u; j < y_size()-1; ++j) {

                if(!obstacle(xpos(i), ypos(j)))
                    (*this)(i, j) = function(i, j);

            }
    }
    else{
        std::cerr << "Obstackle funtion not given" << std::endl;
    }
}

void imn::OGrid::auto_von_Neumann(Condtype type)
{
    auto func = static_cast<vNfunc>(nullptr);

    switch(type){
        case Condtype::PFLOW:
            func = [&,this](unsigned i, unsigned j, Ptype t){ return this->pFlowFunc(i, j, t); };
            break;

        case Condtype::VIFLOW:
            func = [&,this](unsigned i, unsigned j, Ptype t){ return this->viFlowFunc(i, j, t); };
            break;
    }

    if(obstacle){

        for(auto i = 1u; i < x_size()-1; ++i){
            for(auto j = 1u; j < y_size()-1; ++j){

                if( obstacle(xpos(i), ypos(j)) )
                    (*this)(i, j) = func(i, j, con_type(i, j));

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

                if( obstacle(xpos(i), ypos(j)) && con_type(i, j) == pointType)
                    (*this)(i, j) = function(xpos(i), ypos(j));

            }
        }

    }
    else{
        std::cerr << "Obstackle funtion not given" << std::endl;
    }
}

imn::OGrid::Ptype imn::OGrid::con_type(unsigned i, unsigned j) const noexcept
{
    auto left = obstacle(xpos(i-1), ypos(j));
    auto right = obstacle(xpos(i+1), ypos(j));
    auto up = obstacle(xpos(i), ypos(j+1));
    auto down = obstacle(xpos(i), ypos(j-1));

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

double imn::OGrid::pFlowFunc(unsigned i, unsigned j, imn::OGrid::Ptype type) const noexcept
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

double imn::OGrid::viFlowFunc(unsigned i, unsigned j, imn::OGrid::Ptype type) const noexcept
{
    //TODO
    return 0;
}



