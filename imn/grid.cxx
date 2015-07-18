/**
 * @file   grid.hpp
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   July, 2015
 * @brief  Implemenation of the 2D grid class.
 */
 #include "grid.hpp"

 imn::Grid::Grid(int x_min, int x_max, int y_min, int y_max, double dx, double dy, func2d initial, bool widen):
    x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max), dx_(dx), dy_(dy)
 {
    auto size = [](int min, int max, double df)
        {return static_cast<unsigned>(round((max - min) / df));};

    x_size_ = size(x_min_, x_max_, dx_);
    y_size_ = size(y_min_, y_max_, dy_);

     if(widen){
         ++x_size_;
         ++y_size_;
     }

    mtx_ = vMatrix(x_size_ * y_size_, 0);

    if(initial)
        apply_to_all(initial);

 }

void imn::Grid::apply_to_all(func2d function) noexcept
{
    auto i = 0u;
    auto j = 0u;
    for(auto& it : mtx_){
        it = function(xpos(i), ypos(j));

        increm(i, j);
    }
}

void imn::Grid::apply_to_edges(func2d function) noexcept
{
    for(auto i = 0u; i < x_size_; ++i){
        (*this)(i, 0) = function(xpos(i), y_min_);
        (*this)(i, y_size_ -1) = function(xpos(i), y_max_);
    }

    for(auto j = 0u; j < y_size_; ++j){
        (*this)(0, j) = function(x_min_, ypos(j));
        (*this)(x_size_ -1, j) = function(x_max_, ypos(j));
    }
}

void imn::Grid::apply_not_to_edges(func2d function) noexcept
{
    for(auto i = 1u; i < x_size_ -1; ++i)
        for(auto j = 1u; j < y_size_ -1; ++j)
            (*this)(i, j) = function(xpos(i), ypos(j));
}


void imn::Grid::apply_to_single(imn::func2d function, imn::Grid::Edge type) noexcept
{
    switch(type){

        case Edge::left:
            for(auto j = 0u; j < y_size_; ++j)
                (*this)(0, j) = function(x_min_, ypos(j));
            break;

        case Edge::righ:
            for(auto j = 0u; j < y_size_; ++j)
                (*this)(x_size_ -1, j) = function(x_max_, ypos(j));
            break;

        case Edge::up:
            for(auto i = 0u; i < x_size_; ++i)
                (*this)(i, y_size_ -1) = function(xpos(i), y_max_);
            break;

        case Edge::down:
            for(auto i = 0u; i < x_size_; ++i)
                (*this)(i, 0) = function(xpos(i), y_min_);
    }
}

bool imn::Grid::write_to_file(std::ofstream &file, bool zeros) const
{
    if(!file.is_open()) return false;

    auto i = 0u;
    auto j = 0u;
    // I could do that in one for loop, but I think this will be faster
    // In this case, you check zeros only once
    if(zeros){
        for(auto it : mtx_){
            file << xpos(i) << " " << ypos(j) << " " << it << std::endl;

            increm(i, j, &file);
        }
    }
    else{
        // to avoid shitload of blank lines
        auto nn = false;
        for(auto it : mtx_){
            if(it){
                file << xpos(i) << " " << ypos(j) << " " << it << std::endl;
                nn = true;
            }
            if(++j == y_size_){
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
    auto i = 0u;
    auto j = 0u;

    for(auto it : obj.mtx_){
        os << obj.xpos(i) << " " << obj.ypos(j) << " " << it << std::endl;

        obj.increm(i, j, &os);

    }
    return os;
}

void imn::Grid::clear() noexcept
{
    mtx_.assign(x_size_ * y_size_, 0);
}

std::string imn::Grid::grid_to_string() const
{
    std::string s("");
    auto j = 0u;

    for(auto it : mtx_){

        s.append(std::to_string(it) + " ");

        if(!(++j % y_size_))
            s.append("\n");
    }
    return s;
}

inline void imn::Grid::increm(unsigned &i, unsigned &j, std::ostream* stm) const noexcept{
    if(++j == y_size_){
        j = 0;
        ++i;
        if(stm) *stm << std::endl;
    }
}
