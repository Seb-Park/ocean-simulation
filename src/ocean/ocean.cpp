//
// Created by Michael Foiani on 4/9/24.
//

#include <iostream>
#include "ocean.h"

// TODO: make these private instance variables

ocean::ocean()
{
	initial_h = std::vector<std::pair<double, double>>();
    current_h = std::vector<std::pair<double, double>>();

	// initialize the initial height fields
	for (int i = 0; i < N; i++)
	{
        std::pair<double, double> k = index_to_k_vector(i);
        initial_h.push_back(amplitude_0(k));
        current_h.push_back(initial_h[i]);
	}
}

void ocean::updateVertexAmplitudes(double t){
    for (int i = 0; i < N; i++){
        std::pair<double, double> k = index_to_k_vector(i);

        current_h[i] = amplitude_t(t, i, k);
    }
}

/* Maps the 1D k-index into it's 2D waveform vector */
std::pair<double, double> ocean::index_to_k_vector(int index)
{
    // get n' and m' from k_index
    // since k_index goes from 0 <= k_index <= length*width
    int n_prime = index % length;
    int m_prime = index / length;

    double n = (double)n_prime - (double).50*length;
    double m = (double)m_prime - (double).50*width;

    // calculate the k_x and k_z values, according to n and m
    double k_x = (2 * M_PI * n) / L;
    double k_z = (2 * M_PI * m) / L;

    return std::make_pair(k_x, k_z);
}

/* Maps the 1D k-index into it's horizontal x,z position (for passing into shape loader) */
std::pair<double, double> ocean::index_to_horizontal_pos(int index)
{
    // get n' and m' from k_index
    // since k_index goes from 0 <= k_index <= length*width
    int n_prime = index % length;
    int m_prime = index / length;

    double n = (double)n_prime - (double).50*length;
    double m = (double)m_prime - (double).50*width;

    // calculate horizontal pos
    double pos_x = (n * L) / (double) length;
    double pos_z = (m * L) / (double) width;

    return std::make_pair(pos_x, pos_z);
}

/*
 * Randomly generates a complex number by
 * 	importance sampling a Gaussian distribution
 * 	with mean 0 and variance 1.
 *
 * Applies the Box-Muller transform to generate,
 * 	citation: en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
 */
std::pair<double, double> ocean::sample_complex_gaussian
	()
{
	double uniform_1 = (double)rand() / (RAND_MAX);
	double uniform_2 = (double)rand() / (RAND_MAX);

	// set a lower bound on zero to avoid undefined log(0)
	if (uniform_1 == 0)
	{
		uniform_1 = 1e-10;
	}
	if (uniform_2 == 0)
	{
		uniform_2 = 1e-10;
	}

	// real and imaginary parts of the complex number
	double real = sqrt(-2 * log(uniform_1)) * cos(2 * M_PI * uniform_2);
	double imag = sqrt(-2 * log(uniform_1)) * sin(2 * M_PI * uniform_2);

	return std::make_pair(real, imag);
}

/*
 * Generates the Phillips spectrum for a given k-index.
 * See section 4.3 of the paper.
 */
double ocean::phillips_spectrum(std::pair<double, double> k)
{
	// get the k_x, k_z values & amplitude of k
	double k_x = k.first;
	double k_z = k.second;
	double k_magnitude = sqrt(k_x * k_x + k_z * k_z);

	// get the strength of k onto the wind direction
	double k_dot_omega = k_x * omega_wind.first + k_z * omega_wind.second;

	// calculate certain parts of the Phillips formula
	double k_dot_omega_squared = k_dot_omega * k_dot_omega;
	double KL_squared = (k_magnitude * L) * (k_magnitude * L);
	double k_fourth = k_magnitude * k_magnitude * k_magnitude * k_magnitude;

    double l = .1;
	double phillips =
		A // numeric constant
        * exp(-1.0 / KL_squared)
		/ (k_fourth)
        * k_dot_omega_squared
        * exp(-(k_magnitude*k_magnitude*l*l)); 	// consider the small length check as described in the paper

	return phillips;
}

/*
 * Generates the initial (i.e. t = 0) fourier amplitude.
 * See section 4.4 of the paper.
 */
std::pair<double, double> ocean::amplitude_0(std::pair<double, double> k)
{
	std::pair<double, double> xi = sample_complex_gaussian();

	double xi_real = xi.first;
	double xi_imag = xi.second;

    double sqrt_phillips = sqrt(phillips_spectrum(k));

    double inverse_sqrt_2 = 0.707106781187; // this is 1/sqrt(2)
    double real = inverse_sqrt_2 * xi_real * sqrt_phillips;
    double imag = inverse_sqrt_2 * xi_imag * sqrt_phillips;

    return std::make_pair(real, imag);
}

/*
 * Maps the negative k-index from a given k-index,
 * used for the complex conjugate in amplitude_t
 */
int ocean::k_index_to_negative_k_index
	(
		int k_index
	)
{
	int x = k_index % length;
	int z = k_index / length;

	int x_neg = length - x - 1;
	int z_neg = width - z - 1;

	return z_neg * length + x_neg;
}

/*
 * Calculates the dispersion relation for a given k-index.
 * See section 4.2 of the paper.
 */
double ocean::omega_dispersion
	(
		double k_magnitude,
		bool is_shallow
	)
{
	if (is_shallow)
	{
		double tanh_kD = tanh(D * k_magnitude);
        return sqrt(gravity * k_magnitude * tanh_kD);
	}
	else
	{
        return sqrt(gravity * k_magnitude);
	}
}

// where ix is a number made imaginary
std::pair<double, double> ocean::exp_complex(double ix){

    double real = cos(ix);
    double imag = sin(ix);

	return std::make_pair(real, imag);
}

/*
 * Generates the fourier amplitude at time t.
 * See section 4.4 of the paper.
 */
std::pair<double, double> ocean::amplitude_t(double t, int index, std::pair<double, double> k)
{

	// get the initial height (and it's conjugate)
    std::pair<double, double> h_0 = initial_h[index];
    std::pair<double, double> neg_h_0 = amplitude_0(std::make_pair(-k.first, -k.second));
    std::pair<double, double> h_0_conjugate = std::make_pair(h_0.first, h_0.second);;

	// get dispersion from k
	double k_magnitude = sqrt(k.first * k.first + k.second * k.second);
    double omega = omega_dispersion(k_magnitude, M_IS_SHALLOW);

	// calculate the complex exponential terms
	double omega_t = omega * t;
    std::pair<double, double> exp_positive = exp_complex(omega_t);
    std::pair<double, double> exp_negative = exp_complex(-omega_t);

	// add the real and imaginary part together from both h_0 and h_0_conjugate
	double real =
		// h+0 real
        (h_0.first * exp_positive.first)
        + (h_0.second * exp_positive.second)
		// h_0_conjugate real
        + (h_0_conjugate.first * exp_negative.first)
        + (h_0_conjugate.second * exp_negative.second);

	double imag =
		// h_0 imaginary
		h_0.first * exp_positive.second
        + h_0.second * exp_positive.first
		// h_0_conjugate imaginary
		+ h_0_conjugate.first * exp_negative.second
        + h_0_conjugate.second * exp_negative.first;

	return std::make_pair(real, imag);
}

std::vector<Eigen::Vector3f> ocean::get_vertices()
{
	std::vector<Eigen::Vector3f> vertices = std::vector<Eigen::Vector3f>();
	for (int i = 0; i < N; i++)
	{
        std::pair<double, double> k = index_to_horizontal_pos(i);
		double k_x = k.first;
        double k_z = k.second;

        //if (i < length)
        double amplitude = current_h[i].first;
		// double amplitude = sqrt(current_h[i].first * current_h[i].first + current_h[i].second * current_h[i].second);
        // if (i < length) amplitude = initial_h[i].first;


        //std::cout << amplitude << std::endl;

        //std::cout << "k_x: " << k_x << " k_z: " << k_z << " amplitude: " << amplitude << std::endl;

        vertices.emplace_back(k_x*dist_between, amplitude, k_z*dist_between);
	}
	return vertices;
}

std::vector<Eigen::Vector3i> ocean::get_faces()
{
	// connect the vertices into faces
	std::vector<Eigen::Vector3i> faces = std::vector<Eigen::Vector3i>();
	for (int i = 0; i < N; i++)
	{
		int x = i % length;
		int z = i / length;

		// connect the vertices into faces
		if (x < length - 1 && z < width - 1)
		{
			int i1 = i;
			int i2 = i + 1;
			int i3 = i + length;
			int i4 = i + length + 1;

            faces.emplace_back(i2, i1, i3);
            faces.emplace_back(i2, i3, i4);
		}
	}
	return faces;
}
