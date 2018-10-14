#pragma once

namespace conmech
{
namespace stiffness_checker
{

class StiffnessParm
{
 public:
  StiffnessParm():g_(9.80665){}

  ~StiffnessParm(){};

 public:
  // TODO: this radius should be associated with geometry
  /**
   * raidus of the cross section, unit: millimeter
   */
  double radius_;

  /**
   * material density, unit: kg/m^3
   */
  double density_;

  /**
   * unit: MPa
   */
  double youngs_modulus_;

  /**
   * unit: MPa
   */
  double shear_modulus_;

  /**
   * unit: [] (unitless)
   * G = E / 2(1+mu),
   * where G is shear modulus, E is Young's modulus
   * mu is poission ratio
   */
  double poisson_ratio_;

  /**
   * gravity on Earth, unit: m/s^2 = N/kg
   * default value: 9.80665 m/s^2
   * see: https://en.wikipedia.org/wiki/Gravity_of_Earth
   */
  double g_;
};

} // ns stiffness_checker
} // ns conmech