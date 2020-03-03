#pragma once
#include <string>
#include <iostream>

#ifndef NDEBUG
// https://stackoverflow.com/questions/3767869/adding-message-to-assert
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

namespace conmech {

/**
 * @brief numerical zero
 * See: https://github.com/libigl/libigl/blob/master/include/igl/EPS.h
 */
const static double DOUBLE_EPS    = 1.0e-14;
const static double DOUBLE_EPS_SQ = 1.0e-28;
const static double EPSILON = 1.0e-14;

/**
 * @brief unit system used inside conmech. Input unit systems will be converted to these before computation.
 * 
 */
const static std::string LENGTH_UNIT = "meter";
const static std::string ROT_ANGLE_UNIT = "rad";
const static std::string FORCE_UNIT = "kN";

// deduced units
const static std::string MOMENT_UNIT = FORCE_UNIT + "*" + ROT_ANGLE_UNIT;
const static std::string AREA_UNIT         = LENGTH_UNIT + "^2";
const static std::string AREA_INERTIA_UNIT = LENGTH_UNIT + "^4";
const static std::string PRESSURE_UNIT = FORCE_UNIT + "/" + LENGTH_UNIT + "^2";// "kN/m^2" modulus
const static std::string DENSITY_UNIT  = FORCE_UNIT + "/" + LENGTH_UNIT + "^3"; //"kN/m^3"

} // ns conmech