#ifndef POLY_PIC_HELPER_H
#define POLY_PIC_HELPER_H

#include "libWetCloth/Core/MathDefs.h"
#include <Eigen/StdVector>

class PolyPICHelper
{
public:
	PolyPICHelper(const Vector3s& input, const scalar delta_x);
	~PolyPICHelper() = default;

	const scalar Contribution(const int scalar_modes, const VectorXs& coefficients);
	VectorXs CalculateNodeCoefficients(const int scalar_modes, const scalar velocity, const int idx); // Coefficient for the single node, to be summed over all contributing nodes.

private:

	PolyPICHelper() {};

	void CalculateScalarModes();
	void CalculateCoefficientScales(VectorXs& coefficient_scales);

	const scalar G(const scalar input);

	std::vector<scalar> m_scalar_modes;

	scalar m_delta_x_sqr;
	scalar m_inv_delta_x;
	scalar m_inv_delta_x_sqr;

	scalar m_x0;
	scalar m_x1;
	scalar m_x2;

	scalar m_g0;
	scalar m_g1;
	scalar m_g2;
};

#endif