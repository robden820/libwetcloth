#ifndef POLY_PIC_HELPER_H
#define POLY_PIC_HELPER_H

#include "libWetCloth/Core/MathDefs.h"
#include <Eigen/StdVector>

class PolyPICHelper
{
public:
	PolyPICHelper(const Vector3s& node_pos = Vector3s(0), const Vector3s& particle_pos = Vector3s(0), const scalar delta_x = 1.0, const VectorXs& coefficients = VectorXs(0));
	~PolyPICHelper() = default;

	void UpdateValues(const Vector3s& node_pos, const Vector3s& particle_pos, const scalar delta_x, const VectorXs& coefficients = VectorXs(0));

	const scalar Contribution(const int scalar_modes);
	VectorXs CalculateCoefficients(const int scalar_modes, const scalar velocity, const VectorXs& weights, const int idx); // Coefficient for the single node, to be summed over all contributing nodes.

private:

	void CalculateScalarModes();
	void CalculateCoefficientScales();

	const scalar ScalarMode(const int scalar_mode_idx);

	const scalar G(const scalar node_pos, const scalar particle_pos);

	std::vector<scalar> m_scalar_modes;
	std::vector<scalar> m_coefficient_scale;

	scalar m_delta_x_sqr;
	scalar m_inv_delta_x_sqr;

	scalar m_x0;
	scalar m_x1;
	scalar m_x2;

	scalar m_g0;
	scalar m_g1;
	scalar m_g2;

	VectorXs m_coefficients;
};

#endif