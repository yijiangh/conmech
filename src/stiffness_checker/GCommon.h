#pragma once

namespace conmech
{
namespace stiffness_checker
{

#define FILENMAX 128

#ifndef VERSION
#define VERSION "20180921+"
#endif

#ifndef SPT_EPS
#define SPT_EPS  1e-7	// sparse matrix eps
#endif

#ifndef GEO_EPS // geometry eps
#define GEO_EPS  1e-3
#endif

#ifndef STIFF_TOL
#define STIFF_TOL 1.0e-9 // tolerance for stiffness RMS error
#endif

#ifndef MCOND_TOL
#define MCOND_TOL 1.0e12 // tolerance for stiffness matrix condition number
#endif

} // namespace stiffness_checker
} // namespace conmech
