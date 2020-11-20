## Karamba's default rules

- If an element is not assigned with an `Identifier`, then it is automatically assigned with a tag `""`.
- Any `Model` will have a default Material and a default CrossSec. Their corresponding element tags will be `[]` (an empty list of `string`). These are the fall-back properties, i.e. they are used if an element is not assigned with a material/cross sec explicitly by the user.

## References

- [Karamba3d's manual](https://manual.karamba3d.com/). 
- A copy of Karamba3D `1.3.2`'s API is included here for convenience. A more updated version can be obtained [here](https://www.karamba3d.com/support/).