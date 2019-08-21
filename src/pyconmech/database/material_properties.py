"""
Some commonly used material properties.

"""
__all__ = ['MATERIAL_PROPERTY_DATABASE']

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
        "radius": 2.5,
        "radius_unit": "centimeter"
    }, 
    'PLA' : {
        "_source": "https://www.makeitfrom.com/material-properties/Polylactic-Acid-PLA-Polylactide",
        "material_name": "PLA",
        "youngs_modulus": 350,
        "youngs_modulus_unit": "kN/cm2",
        "shear_modulus": 240,
        "shear_modulus_unit": "kN/cm2",
        "density": 0.000012748,
        "density_unit": "kN/m3",
        "poisson_ratio": -0.27083,
        "radius": 2.5,
        "radius_unit": "centimeter"
    }, 
}