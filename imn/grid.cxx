/**
 * \file   grid.hpp
 * \Author Bartłomiej Meder (bartem93@gmail.com)
 * \date   July, 2015
 * \brief  Implemenation of 2D grid class.
 */
 #include "grid.hpp"

 imn::Grid::Grid(int x_min, int x_max, int y_min, int y_max, double dx, double dy, func2d initial):
    _x_min(x_min), _x_max(x_max), _y_min(y_min), _y_max(y_max), _dx(dx), _dy(dy)
 {
    auto size = [](int min, int max, double df)
        {return static_cast<std::size_t>(round((max - min) / df));};

    _x_size = size(_x_min, _x_max, _dx);
    _y_size = size(_y_min, _y_max, _dy);

    _mtx = vMatrix(_x_size * _y_size, 0);

    if(initial)
        apply_to_all(initial);

 }

void imn::Grid::apply_to_all(func2d function) noexcept
{
    auto i = 0;
    auto j = 0;
    for(auto& it : _mtx){
       it = function(_x_min + i*_dx, _y_min + j*_dy);
       if(++j == _y_size){
           j = 0;
           ++i;
       }
    }
}

void imn::Grid::apply_to_edge(func2d function) noexcept
{
    for(std::size_t i = 0; i < _x_size; ++i){
        (*this)(i, 0) = function(_x_min + i*_dx, _y_min);
        (*this)(i, _y_size-1) = function(_x_min + i*_dx, _y_max);
    }

    for(std::size_t j = 0; j < _y_size; ++j){
        (*this)(0, j) = function(_x_min, _y_min + j*_dy);
        (*this)(_x_size-1, j) = function(_x_max, _y_min + j*_dy);
    }
}

void imn::Grid::apply_not_to_edge(func2d function) noexcept
{
    for(std::size_t i = 1; i < _x_size-1; ++i)
        for(std::size_t j = 1; j < _y_size-1; ++j)
            (*this)(i, j) = function(_x_min + i*_dx, _y_min + j*_dy);
}

bool imn::Grid::write_to_file(std::ofstream &file, bool zeros) const
{
    if(!file.is_open()) return false;

    auto i = 0;
    auto j = 0;
    // I could do that in one for loop, but I think this will be faster
    // In this case, you check zeros only once
    if(zeros){
        for(auto it : _mtx){
            file << _x_min + i*_dx << " " << _y_min + j*_dy << " " << it << std::endl;
            if(++j == _y_size){
                j = 0;
                ++i;
                file << std::endl;
            }
        }
    }
    else{
        // to avoid shitload of blank lines
        auto nn = false;
        for(auto it : _mtx){
            if(it){
                file << _x_min + i*_dx << " " << _y_min + j*_dy << " " << it << std::endl;
                nn = true;
            }
            if(++j == _y_size){
                j = 0;
                ++i;
                if(nn){ file << std::endl; nn = false; }
            }
        }
    }

    return true;
}

std::ostream& imn::operator<<(std::ostream& os, const imn::Grid& obj)
{
    auto i = 0;
    auto j = 0;

    for(auto it : obj._mtx){
        os << obj._x_min + i*obj._dx << " " << obj._y_min + j*obj._dy << " " << it << std::endl;

        if(++j == obj._y_size){
            j = 0;
            ++i;
            os << std::endl;
        }

    }
    return os;
}

void imn::Grid::clear() noexcept{
    for(auto& it : _mtx)
        it = 0;
}

std::string imn::Grid::grid_to_string() const {
    std::string s("");
    auto j=0;
    for(auto it: _mtx){
        s.append(std::to_string(it) + " ");
        if(++j == _y_size){
            j = 0;
            s.append("\n");
        }
    }
    return s;
}
