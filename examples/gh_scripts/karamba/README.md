# Karamba Problem Instance Generator

This is a [Grasshopper](https://www.grasshopper3d.com/) script on the [Rhinoceros](https://www.rhino3d.com/) platform, 
used for creating problem instances for `conmech`.

The Grasshopper script depends on [Karamba3d](https://www.karamba3d.com/) for assembling the model (elements, fixities, loads, materials, cross sections, etc.). 
The free version of Karamba3d would suffice since we don't use its analysis components.

We uses the `karamba_io` modules from `pyconmech` to do two-way conversions between the Karamba models and `json` files. The `karamba_io.py` and `io_base.py` has minimal dependencies and can be used separately from `conmech`.

<img src="./images/karamba_model_export.png">

## Karamba's default rules

- If an element is not assigned with an `Identifier`, then it is automatically assigned with a tag `""`.
- Any `Model` will have a default Material and a default CrossSec. Their corresponding element tags will be `[]` (an empty list of `string`). These are the fall-back properties, i.e. they are used if an element is not assigned with a material/cross sec explicitly by the user.

## References

- [Karamba3d's manual](https://manual.karamba3d.com/). 
- A copy of Karamba3D `1.3.2`'s API is included here for convenience. A more updated version can be obtained [here](https://www.karamba3d.com/support/).