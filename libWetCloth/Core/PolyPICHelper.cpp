#include "PolyPICHelper.h"
#include "MathUtilities.h"

PolyPICHelper::PolyPICHelper(const Vector3s& input, const scalar delta_x)
	: m_delta_x_sqr(delta_x * delta_x), m_inv_delta_x(1 / delta_x), m_inv_delta_x_sqr(m_inv_delta_x * m_inv_delta_x),
	m_x0(input.x()), m_x1(input.y()), m_x2(input.z()),
	m_g0(G(m_x0)), m_g1(G(m_x1)), m_g2(G(m_x2))
{
	m_scalar_modes.resize(27);

	CalculateScalarModes();
}

void PolyPICHelper::CalculateScalarModes()
{
	m_scalar_modes[0] = 1;
	m_scalar_modes[1] = m_x0;
	m_scalar_modes[2] = m_x1;
	m_scalar_modes[3] = m_x2;
	m_scalar_modes[4] = m_x0 * m_x1;
	m_scalar_modes[5] = m_x0 * m_x2;
	m_scalar_modes[6] = m_x1 * m_x2;
	m_scalar_modes[7] = m_x0 * m_x1 * m_x2;
	m_scalar_modes[8] = m_g0;
	m_scalar_modes[9] = m_g1;
	m_scalar_modes[10] = m_g2;
	m_scalar_modes[11] = m_g0 * m_g1;
	m_scalar_modes[12] = m_g1 * m_g2;
	m_scalar_modes[13] = m_g0 * m_g2;
	m_scalar_modes[14] = m_g0 * m_g1 * m_g2;
	m_scalar_modes[15] = m_g0 * m_x1;
	m_scalar_modes[16] = m_g0 * m_x2;
	m_scalar_modes[17] = m_g0 * m_x1 * m_x2;
	m_scalar_modes[18] = m_g1 * m_x0;
	m_scalar_modes[19] = m_g1 * m_x2;
	m_scalar_modes[20] = m_g1 * m_x0 * m_x2;
	m_scalar_modes[21] = m_g2 * m_x0;
	m_scalar_modes[22] = m_g2 * m_x1;
	m_scalar_modes[23] = m_g2 * m_x0 * m_x1;
	m_scalar_modes[24] = m_g0 * m_g1 * m_x2;
	m_scalar_modes[25] = m_g0 * m_g2 * m_x1;
	m_scalar_modes[26] = m_g1 * m_g2 * m_x0;
}

void PolyPICHelper::CalculateCoefficientScales(VectorXs& coefficient_scales)
{
	assert(coefficient_scales.size() >= 27);

	// See tech_doc for following constants.
	// Constants are calculated as 1 / presented_value.

	const scalar inv_delta_x_pow_4 = m_inv_delta_x_sqr * m_inv_delta_x_sqr;
	const scalar inv_delta_x_pow_6 = inv_delta_x_pow_4 * m_inv_delta_x_sqr;

	coefficient_scales[0] = 1;
	coefficient_scales[1] = m_inv_delta_x_sqr * 4;
	coefficient_scales[2] = m_inv_delta_x_sqr * 4;
	coefficient_scales[3] = m_inv_delta_x_sqr * 4;
	coefficient_scales[4] = inv_delta_x_pow_4 * 16;
	coefficient_scales[5] = inv_delta_x_pow_4 * 16;
	coefficient_scales[6] = inv_delta_x_pow_4 * 16;
	coefficient_scales[7] = inv_delta_x_pow_6 * 64;
	coefficient_scales[8] = m_inv_delta_x_sqr * 0.5;
	coefficient_scales[9] = m_inv_delta_x_sqr * 0.5;
	coefficient_scales[10] = m_inv_delta_x_sqr * 0.5;
	coefficient_scales[11] = inv_delta_x_pow_4 * 0.25;
	coefficient_scales[12] = inv_delta_x_pow_4 * 0.25;
	coefficient_scales[13] = inv_delta_x_pow_4 * 0.25;
	coefficient_scales[14] = inv_delta_x_pow_6 * 0.125;
	coefficient_scales[15] = inv_delta_x_pow_4 * 2;
	coefficient_scales[16] = inv_delta_x_pow_4 * 2;
	coefficient_scales[17] = inv_delta_x_pow_6 * 8;
	coefficient_scales[18] = inv_delta_x_pow_4 * 2;
	coefficient_scales[19] = inv_delta_x_pow_4 * 2;
	coefficient_scales[20] = inv_delta_x_pow_6 * 8;
	coefficient_scales[21] = inv_delta_x_pow_4 * 2;
	coefficient_scales[22] = inv_delta_x_pow_4 * 2;
	coefficient_scales[23] = inv_delta_x_pow_6 * 8;
	coefficient_scales[24] = inv_delta_x_pow_6;
	coefficient_scales[25] = inv_delta_x_pow_6;
	coefficient_scales[26] = inv_delta_x_pow_6;
}

const scalar PolyPICHelper::Contribution(const int scalar_modes, const VectorXs& coefficients)
{
	assert(scalar_modes > 0 && scalar_modes <= 27);
	assert(coefficients.size() >= scalar_modes);

	scalar mode_sum = 0.0;

	for (int mode = 0; mode < scalar_modes; mode++)
	{
		mode_sum += m_scalar_modes[mode] * coefficients[mode];
	}

	return mode_sum;
}

VectorXs PolyPICHelper::CalculateNodeCoefficients(const int scalar_modes, const scalar velocity, const int idx)
{
	assert(scalar_modes > 0 && scalar_modes <= 27);
	assert(idx >= 0 && idx < 27);

	VectorXs coefficients(27);
	VectorXs coefficient_scales(27);

	CalculateCoefficientScales(coefficient_scales);

	const scalar weight_0 = mathutils::quad_kernel(m_x0 * m_inv_delta_x); // weight for x axis.
	const scalar weight_1 = mathutils::quad_kernel(m_x1 * m_inv_delta_x); // y axis.
	const scalar weight_2 = mathutils::quad_kernel(m_x2 * m_inv_delta_x); // z axis.
	const scalar weight = weight_0 * weight_1 * weight_2;

	const int idx_0 = std::fmod(idx, 3.0);
	const int idx_1 = std::fmod(std::floor(idx / 3.0), 3.0);
	const int idx_2 = std::fmod(std::floor(idx / 9.0), 3.0);

	const scalar exp_0 = std::pow(-2.0, idx_0 % 2);
	const scalar exp_1 = std::pow(-2.0, idx_1 % 2);
	const scalar exp_2 = std::pow(-2.0, idx_2 % 2);

	coefficients[0] = weight;
	coefficients[1] = weight * m_x0;
	coefficients[2] = weight * m_x1;
	coefficients[3] = weight * m_x2;
	coefficients[4] = weight * m_x0 * m_x1;
	coefficients[5] = weight * m_x1 * m_x2;
	coefficients[6] = weight * m_x0 * m_x2;
	coefficients[7] = weight * m_x0 * m_x1 * m_x2;
	coefficients[8] = weight_1 * weight_2 * exp_0;
	coefficients[9] = weight_0 * weight_2 * exp_1;
	coefficients[10] = weight_0 * weight_1 * exp_2;
	coefficients[11] = weight_2 * exp_0 * exp_1;
	coefficients[12] = weight_1 * exp_0 * exp_2;
	coefficients[13] = weight_0 * exp_1 * exp_2;
	coefficients[14] = exp_0 * exp_1 * exp_2;
	coefficients[15] = weight_1 * weight_2 * m_x1 * exp_0;
	coefficients[16] = weight_1 * weight_2 * m_x2 * exp_0;
	coefficients[17] = weight_1 * weight_2 * m_x1 * m_x2 * exp_0;
	coefficients[18] = weight_0 * weight_2 * m_x0 * exp_1;
	coefficients[19] = weight_0 * weight_2 * m_x2 * exp_1;
	coefficients[20] = weight_0 * weight_2 * m_x0 * m_x2 * exp_1;
	coefficients[21] = weight_0 * weight_1 * m_x0 * exp_2;
	coefficients[22] = weight_0 * weight_1 * m_x1 * exp_2;
	coefficients[23] = weight_0 * weight_1 * m_x0 * m_x1 * exp_2;
	coefficients[24] = weight_2 * m_x2 * exp_0 * exp_1;
	coefficients[25] = weight_1 * m_x1 * exp_0 * exp_2;
	coefficients[26] = weight_0 * m_x0 * exp_1 * exp_2;
	
	for (int i = 0; i < scalar_modes; i++)
	{
		coefficients[i] *= coefficient_scales[i] * velocity;
	}
	
	return coefficients.segment(0, scalar_modes);
}

const scalar PolyPICHelper::G(const scalar input)
{
	const scalar weight = mathutils::quad_kernel(input * m_inv_delta_x);

	const scalar a = weight * weight;
	const scalar b = input * (m_delta_x_sqr - 4.0 * input * input) * weight * m_inv_delta_x_sqr;
	const scalar c = m_delta_x_sqr * 0.25;

	return a - b - c;
}
