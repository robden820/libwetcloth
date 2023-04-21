#include "PolyPICHelper.h"
#include "MathUtilities.h"

PolyPICHelper::PolyPICHelper(const Vector3s& node_pos, const Vector3s& particle_pos, const scalar delta_x)
	: m_inv_delta_x(1.0 / delta_x), m_delta_x_sqr(delta_x * delta_x), m_inv_delta_x_sqr(m_inv_delta_x * m_inv_delta_x),
	m_diff(node_pos - particle_pos),
	m_x0(m_diff.x()), m_x1(m_diff.y()), m_x2(m_diff.z()),
	m_g0(G(node_pos.x(),particle_pos.x())), m_g1(G(node_pos.y(),particle_pos.y())), m_g2(G(node_pos.z(),particle_pos.z()))
{
	m_scalar_modes.resize(27);
	m_coefficient_scale.resize(27);

	CalculateScalarModes();
	CalculateCoefficientScales();
}

void PolyPICHelper::UpdateValues(const Vector3s& node_pos, const Vector3s& particle_pos, const scalar delta_x)
{
	assert(delta_x != 0.0);

	if (delta_x * delta_x != m_delta_x_sqr)
	{
		m_inv_delta_x = 1.0 / delta_x;
		m_delta_x_sqr = delta_x * delta_x;
		m_inv_delta_x_sqr = m_inv_delta_x * m_inv_delta_x;

		CalculateCoefficientScales();
	}

	m_diff = (node_pos - particle_pos);

	m_x0 = m_diff.x();
	m_x1 = m_diff.y();
	m_x2 = m_diff.z();

	m_g0 = G(m_x0, 0);
	m_g1 = G(m_x1, 1);
	m_g2 = G(m_x2, 2);

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

void PolyPICHelper::CalculateCoefficientScales()
{
	const scalar inv_delta_x_pow_4 = m_inv_delta_x_sqr * m_inv_delta_x_sqr;
	const scalar inv_delta_x_pow_6 = inv_delta_x_pow_4 * m_inv_delta_x_sqr;

	m_coefficient_scale[0] = 1;
	m_coefficient_scale[1] = m_inv_delta_x_sqr * 4;
	m_coefficient_scale[2] = m_inv_delta_x_sqr * 4;
	m_coefficient_scale[3] = m_inv_delta_x_sqr * 4;
	m_coefficient_scale[4] = inv_delta_x_pow_4 * 16;
	m_coefficient_scale[5] = inv_delta_x_pow_4 * 16;
	m_coefficient_scale[6] = inv_delta_x_pow_4 * 16;
	m_coefficient_scale[7] = inv_delta_x_pow_6 * 64;
	m_coefficient_scale[8] = m_inv_delta_x_sqr * 0.5;
	m_coefficient_scale[9] = m_inv_delta_x_sqr * 0.5;
	m_coefficient_scale[10] = m_inv_delta_x_sqr * 0.5;
	m_coefficient_scale[11] = inv_delta_x_pow_4 * 0.25;
	m_coefficient_scale[12] = inv_delta_x_pow_4 * 0.25;
	m_coefficient_scale[13] = inv_delta_x_pow_4 * 0.25;
	m_coefficient_scale[14] = inv_delta_x_pow_6 * 0.125;
	m_coefficient_scale[15] = inv_delta_x_pow_4 * 2;
	m_coefficient_scale[16] = inv_delta_x_pow_4 * 2;
	m_coefficient_scale[17] = inv_delta_x_pow_6 * 8;
	m_coefficient_scale[18] = inv_delta_x_pow_4 * 2;
	m_coefficient_scale[19] = inv_delta_x_pow_4 * 2;
	m_coefficient_scale[20] = inv_delta_x_pow_6 * 8;
	m_coefficient_scale[21] = inv_delta_x_pow_4 * 2;
	m_coefficient_scale[22] = inv_delta_x_pow_4 * 2;
	m_coefficient_scale[23] = inv_delta_x_pow_6 * 8;
	m_coefficient_scale[24] = inv_delta_x_pow_6;
	m_coefficient_scale[25] = inv_delta_x_pow_6;
	m_coefficient_scale[26] = inv_delta_x_pow_6;
}

const scalar PolyPICHelper::Contribution(const int scalar_modes, const VectorXs& coefficients)
{
	assert(scalar_modes > 0 && scalar_modes < 27);
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
	VectorXs coefficients(27);

	const scalar weight_0 = mathutils::quad_kernel(m_x0 * m_inv_delta_x); // weight for x axis interpolation.
	const scalar weight_1 = mathutils::quad_kernel(m_x1 * m_inv_delta_x); // y axis.
	const scalar weight_2 = mathutils::quad_kernel(m_x2 * m_inv_delta_x); // z axis.
	const scalar weight = weight_0 * weight_1 * weight_2;

	const int exp_0 = std::fmod(idx, 3.0);
	const int exp_1 = std::fmod(std::floor(idx / 3.0), 3.0);
	const int exp_2 = std::fmod(std::floor(idx / 9.0), 3.0);

	const scalar two_pow_0 = std::pow(-2, exp_0 - 1 % 2);
	const scalar two_pow_1 = std::pow(-2, exp_1 - 1 % 2);
	const scalar two_pow_2 = std::pow(-2, exp_2 - 1 % 2);

	coefficients[0] = weight;
	coefficients[1] = weight * m_x0;
	coefficients[2] = weight * m_x1;
	coefficients[3] = weight * m_x2;
	coefficients[4] = weight * m_x0 * m_x1;
	coefficients[5] = weight * m_x1 * m_x2;
	coefficients[6] = weight * m_x0 * m_x2;
	coefficients[7] = weight * m_x0 * m_x1 * m_x2;
	coefficients[8] = weight_1 * weight_2 * two_pow_0;
	coefficients[9] = weight_0 * weight_2 * two_pow_1;
	coefficients[10] = weight_0 * weight_1 * two_pow_2;
	coefficients[11] = weight_2 * two_pow_0 * two_pow_1;
	coefficients[12] = weight_1 * two_pow_0 * two_pow_2;
	coefficients[13] = weight_0 * two_pow_1 * two_pow_2;
	coefficients[14] = two_pow_0 * two_pow_1 * two_pow_2;
	coefficients[15] = weight_1 * weight_2 * m_x1 * two_pow_0;
	coefficients[16] = weight_1 * weight_2 * m_x2 * two_pow_0;
	coefficients[17] = weight_1 * weight_2 * m_x1 * m_x2 * two_pow_0;
	coefficients[18] = weight_0 * weight_2 * m_x0 * two_pow_1;
	coefficients[19] = weight_0 * weight_2 * m_x2 * two_pow_1;
	coefficients[20] = weight_0 * weight_2 * m_x0 * m_x2 * two_pow_1;
	coefficients[21] = weight_0 * weight_1 * m_x0 * two_pow_2;
	coefficients[22] = weight_0 * weight_1 * m_x1 * two_pow_2;
	coefficients[23] = weight_0 * weight_1 * m_x0 * m_x1 * two_pow_2;
	coefficients[24] = weight_2 * m_x2 * two_pow_0 * two_pow_1;
	coefficients[25] = weight_1 * m_x1 * two_pow_0 * two_pow_2;
	coefficients[26] = weight_0 * m_x0 * two_pow_1 * two_pow_2;
	
	for (unsigned int i = 0; i < 27; i++)
	{
		coefficients[i] *= m_coefficient_scale[i] * velocity;
	}
	
	return coefficients.segment(0, scalar_modes);
}

const scalar PolyPICHelper::ScalarMode(const int scalar_mode_idx)
{
	assert(scalar_mode_idx >= 0 && scalar_mode_idx < 27);

	return m_scalar_modes[scalar_mode_idx];
}

const scalar PolyPICHelper::G(const scalar node_pos, const scalar particle_pos)
{
	const scalar diff = node_pos - particle_pos;

	const scalar a = diff * diff;
	const scalar b = (particle_pos * (m_delta_x_sqr - 4 * particle_pos * particle_pos) * diff) * m_inv_delta_x_sqr;
	const scalar c = m_delta_x_sqr * 0.25;

	return a - b - c;
}
