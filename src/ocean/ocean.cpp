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
		initial_h.push_back(amplitude_0(i));
        current_h.push_back(initial_h[i]);
	}
}

void ocean::updateVertexAmplitudes(double t){
    for (int i = 0; i < N; i++){
        current_h[i] = amplitude_t(t, i);
    }
}

/* Maps the 1D k-index into it's 2D waveform vector */
std::pair<double, double> ocean::k_index_to_k_vector(int k_index)
{
	// get the x and z indices
	int x = k_index % length;
	int z = k_index / length;

	// calculate the k_x and k_z values, according to the length
	double k_x = (2 * M_PI * x - N) / (double) length;
	double k_z = (2 * M_PI * z - N) / (double) width;

	return std::make_pair(k_x, k_z);
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
double ocean::phillips_spectrum
	(
	int k_index
	)
{
	// get the k_x, k_z values & amplitude of k
	std::pair<double, double> k = k_index_to_k_vector(k_index);
	double k_x = k.first;
	double k_z = k.second;
	double k_magnitude = sqrt(k_x * k_x + k_z * k_z);

	// calculate L, the max wave size from wind speed V
	double L = (V * V) / 9.81;
	// get the strength of k onto the wind direction
	double k_dot_omega = k_x * omega_wind.first + k_z * omega_wind.second;

	// calculate certain parts of the Phillips formula
	double k_dot_omega_squared = k_dot_omega * k_dot_omega;
	double KL_squared = (k_magnitude * L) * (k_magnitude * L);
	double k_fourth = k_magnitude * k_magnitude * k_magnitude * k_magnitude;

	double phillips =
		A // numeric constant
		* exp(-1 / KL_squared)
		/ (k_fourth)
		* k_dot_omega_squared;

	// TODO: consider the small length check as described in the paper
	return phillips;
}

/*
 * Generates the initial (i.e. t = 0) fourier amplitude.
 * See section 4.4 of the paper.
 */
std::pair<double, double> ocean::amplitude_0
	(
		int k_index
	)
{
	std::pair<double, double> xi = sample_complex_gaussian();

	double xi_real = xi.first;
	double xi_imag = xi.second;

	double sqrt_phillips = sqrt(phillips_spectrum(k_index));

	double real = (1 / sqrt(2)) * xi_real * sqrt_phillips;
	double imag = (1 / sqrt(2)) * xi_imag * sqrt_phillips;

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
		return sqrt(9.81 * k_magnitude * tanh_kD);
	}
	else
	{
		return sqrt(9.81 * k_magnitude);
	}
}

std::pair<double, double> ocean::exp_complex
	(
		std::pair<double, double> z
	)
{
	double real = exp(z.first) * cos(z.second);
	double imag = exp(z.first) * sin(z.second);

	return std::make_pair(real, imag);
}

/*
 * Generates the fourier amplitude at time t.
 * See section 4.4 of the paper.
 */
std::pair<double, double> ocean::amplitude_t
	(
		double t,
		int k_index
	)
{
	// get the initial height (and it's conjugate)
	std::pair<double, double> h_0 = initial_h[k_index];
	std::pair<double, double> h_0_conjugate = initial_h[k_index_to_negative_k_index(k_index)];
	h_0_conjugate = std::make_pair(h_0_conjugate.first, -h_0_conjugate.second);

	// get dispersion from k
	std::pair<double, double> k = k_index_to_k_vector(k_index);
	double k_magnitude = sqrt(k.first * k.first + k.second * k.second);
	double omega = omega_dispersion(k_magnitude, true);

	// calculate the complex exponential terms
	double omega_t = omega * t;
	std::pair<double, double> exp_positive = exp_complex(std::make_pair(0, omega_t));
	std::pair<double, double> exp_negative = exp_complex(std::make_pair(0, -omega_t));

	// add the real and imaginary part together from both h_0 and h_0_conjugate
	double real =
		// h+0 real
        (h_0.first * exp_positive.first)
        - (h_0.second * exp_positive.second)
		// h_0_conjugate real
        + (h_0_conjugate.first * exp_negative.first)
        - (h_0_conjugate.second * exp_negative.second);

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
		std::pair<double, double> k = k_index_to_k_vector(i);
		double k_x = k.first;
		double k_z = k.second;

        //if (i < length)
        double amplitude = current_h[i].first;
		// double amplitude = sqrt(current_h[i].first * current_h[i].first + current_h[i].second * current_h[i].second);
        // if (i < length) amplitude = initial_h[i].first;


        //if (i==2) std::cout << amplitude << std::endl;

        //std::cout << "k_x: " << k_x << " k_z: " << k_z << " amplitude: " << amplitude << std::endl;

		vertices.emplace_back(k_x, amplitude, k_z);
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

//
//std::vector<Eigen::Vector2d> ocean_alt::fast_fft
//	(
//		std::vector<Eigen::Vector2d> h
//	)
//{
//	int N = h.size();
//	std::vector<Eigen::Vector2d> H = std::vector<Eigen::Vector2d>(N);
//
//	if (N == 1)
//	{
//		H[0] = h[0];
//		return H;
//	}
//	else
//	{
//		std::vector<Eigen::Vector2d> even = std::vector<Eigen::Vector2d>(N / 2);
//		std::vector<Eigen::Vector2d> odd = std::vector<Eigen::Vector2d>(N / 2);
//
//		for (int i = 0; i < N / 2; i++)
//		{
//			even[i] = h[2 * i];
//			odd[i] = h[2 * i + 1];
//		}
//
//		std::vector<Eigen::Vector2d> even_fft = fast_fft(even);
//		std::vector<Eigen::Vector2d> odd_fft = fast_fft(odd);
//
//		for (int i = 0; i < N / 2; i++)
//		{
//			Eigen::Vector2d x_vector = m_waveIndexConstants[i].base_horiz_pos;
//			Eigen::Vector2d k_vector = m_waveIndexConstants[i].k_vector;
//
//			Eigen::Vector2d h_tilda_prime = h[i]; // vector(real, imag)
//
//			double k_dot_x = x_vector.dot(k_vector);
//			Eigen::Vector2d omega = complex_exp(-k_dot_x); // vector(real, imag)
//			Eigen::Vector2d omega_times_odd =
//				{
//					omega[0] * odd_fft[i][0] - omega[1] * odd_fft[i][1],
//					omega[0] * odd_fft[i][1] + omega[1] * odd_fft[i][0]
//				};
//
//			H[i] = Eigen::Vector2d
//				(
//					even_fft[i][0] + omega_times_odd[0],
//					even_fft[i][1] + omega_times_odd[1]
//				);
//			H[i + N / 2] = Eigen::Vector2d
//				(
//					even_fft[i][0] - omega_times_odd[0],
//					even_fft[i][1] - omega_times_odd[1]
//				);
//		}
//
//		return H;
//	}
//}
