#include "complex.h"


// squre root of the complex number z 
complex  sqrt(const complex& z)
{
	if (z == 0.0) return 0;
	double	r = fabs(real(z)), i = fabs(imag(z)), t,
		w = (r >= i) ? (t = i / r, sqrt(r) * sqrt(0.5 * (1.0 + sqrt(1.0 + t * t)))) :
		(t = r / i, sqrt(i) * sqrt(0.5 * (1.0 + sqrt(1.0 + t * t))));

	return real(z) > 0 ? complex(w, imag(z) / (2.0 * w)) :
		((t = imag(z) >= 0.0 ? w : -w), complex(imag(z) / (2.0 * t), t));
}
//------------------------------------------------------------------------------

complex  exp_im(double x)
{
	return complex(cos(x), sin(x));
};
//------------------------------------------------------------------------------
complex zeroComplex() { return complex(0.0, 0.0); }