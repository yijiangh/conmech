"""
Some commonly used material properties.

"""
__all__ = ['MATERIAL_PROPERTY_DATABASE']

"""material properties library: Steel, PLA
"""
MATERIAL_PROPERTY_DATABASE = {
    'Steel-S235' : {
        "material_name": "Steel-S235",
        "youngs_modulus": 21000.0,
        "youngs_modulus_unit": "kN/cm2",
        "shear_modulus": 8076.0,
        "shear_modulus_unit": "kN/cm2",
        "density": 78.5,
        "density_unit": "kN/m3",
        "poisson_ratio": 0.300149,
        "cross_sec_area": 19.634954084936208,
        "cross_sec_area_unit": "centimeter^2",
        "_cross_sec_area_note": "round_shape, radius 2.5 cm",
        "Jx": 61.35923151542565,
        "Jx_unit": "centimeter^4",
        "_Jx_note": "0.5 * pi * r^4",
        "Iy": 30.679615757712824,
        "Iy_unit": "centimeter^4",
        "_Iy_note": "0.25 * pi * r^4",
        "Iz": 30.679615757712824,
        "Iz_unit": "centimeter^4",
        "_Iz_note": "0.25 * pi * r^4"
     }, 
    'PLA' : {
        "_source": "https://www.makeitfrom.com/material-properties/Polylactic-Acid-PLA-Polylactide",
        "material_name": "PLA",
        "youngs_modulus": 350,
        "youngs_modulus_unit": "kN/cm2",
        "shear_modulus": 240,
        "shear_modulus_unit": "kN/cm2",
        "tensile_yeild_stress": 3.6,
        "tensile_yeild_stress_unit": "kN/cm2",
        "density": 12.2582,
        "density_unit": "kN/m3",
        "poisson_ratio": -0.27083,
        "cross_sec_area": 0.07068583470577035,
        "cross_sec_area_unit": "centimeter^2",
        "_cross_sec_area_note": "round_shape, radius 0.15 cm",
        "Jx": 0.0007952156404399163,
        "Jx_unit": "centimeter^4",
        "_Jx_note": "0.5 * pi * r^4",
        "Iy": 0.00039760782021995816,
        "Iy_unit": "centimeter^4",
        "_Iy_note": "0.25 * pi * r^4",
        "Iz": 0.00039760782021995816,
        "Iz_unit": "centimeter^4",
        "_Iz_note": "0.25 * pi * r^4"
  }, 
}