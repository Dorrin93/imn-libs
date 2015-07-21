/**
 * @file   p_flow.hpp
 * @Author Bartłomiej Meder (bartem93@gmail.com)
 * @date   July, 2015
 * @brief  TODO
 *
 * TODO
 */
#ifndef __IMN_P_FLOW__
#define __IMN_P_FLOW__
#include "ogrid.hpp"
#include "utils.hpp"
#include <fstream>

/**
 * @addtogroup imn
 * @{
 */

//* Engineering Numerical Methods functions
namespace imn {

    void stream_lines(OGrid& grd, double u0, std::ofstream& file, double eps = 1e-6, bool w_grid = true);

}

/** @} End of group */

#endif
