#pragma once

#define FILENMAX 128

#ifndef VERSION
#define VERSION "20160826+"
#endif

#ifndef F_PI
#define F_PI 3.14159265358979323846264338327950288419716939937510
#endif

#ifndef SPT_EPS
#define SPT_EPS  0.0000001	// sparse matrix eps
#endif

#ifndef  GEO_EPS			// geometry eps
#define GEO_EPS  0.001
#endif

#ifndef STIFF_TOL
#define STIFF_TOL 1.0e-9	// tolerance for stiffness RMS error
#endif

#ifndef MCOND_TOL
#define MCOND_TOL 1.0e12	// tolerance for stiffness matrix condition number
#endif

// Zvert=1: Z axis is vertical... rotate about Y-axis, then rotate about Z-axis
// Zvert=0: Y axis is vertical... rotate about Z-axis, then rotate about Y-axis
#define Zvert 1	

#endif /* FIBERPRINT_COMMON_H */
