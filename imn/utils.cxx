/**
 * @file   utils.cxx
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   July, 2015
 * @brief  Useful funtions implementations.
 */
#include "utils.hpp"

std::string imn::double_to_string(double x)
{
    auto s_x = std::to_string(x);

    s_x.erase(s_x.find_last_not_of('0') + 1, std::string::npos);

    if(s_x[s_x.length()-1] == '.')
        s_x.append("0");

    return s_x;
}