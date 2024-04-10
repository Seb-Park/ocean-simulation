//
// Created by Michael Foiani on 4/9/24.
//

#ifndef OCEAN_H
#define OCEAN_H

#include <vector>
#include <utility>
#include <Eigen/Dense>


class ocean
{
public:
	ocean();
    void updateVertexAmplitudes(double t);
    std::vector<Eigen::Vector3f> get_vertices();
    std::vector<Eigen::Vector3i> get_faces();




private:

    const int length = 256; // length of grid
    const int width = 256; // width of grid
	const int N = length * width; // total number of grid points

    const double A = 30; // numeric constant for the Phillips spectrum
    const double V = 0.5; // wind speed
	const std::pair<double, double> omega_wind
		= std::make_pair(1.0, 0.0); // wind direction

	std::vector<std::pair<double, double>> initial_h; // initial height fields for each K
    std::vector<std::pair<double, double>> current_h; // current height fields for each K


	const double D = 1; // Depth below mean water level (for dispersion relation)


	std::pair<double, double> k_index_to_k_vector
	(
		int k_index
	);
	std::pair<double, double> sample_complex_gaussian
		();
	double phillips_spectrum
	(
		int k_index
	);
	std::pair<double, double> amplitude_0
	(
		int k_index
	);
	int k_index_to_negative_k_index
	(
		int k_index
	);
	double omega_dispersion
	(
		double k_magnitude,
		bool is_shallow=false
	);
	std::pair<double, double> exp_complex
	(
		std::pair<double, double> z
	);
	std::pair<double, double> amplitude_t
	(
		double t,
		int k_index
	);
};


#endif //OCEAN_H
