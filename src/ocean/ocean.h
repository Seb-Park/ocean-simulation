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

    const bool M_IS_SHALLOW = false;

    const int length = 512; // length of grid
    const int width = 512; // width of grid
	const int N = length * width; // total number of grid points
    const double dist_between = 10.0; // distance between grid points

    const double A = 30.0; // numeric constant for the Phillips spectrum
    const double V = 100; // wind speed
    const double gravity = 9.81;
    const double L = V*V/gravity;
	const std::pair<double, double> omega_wind
		= std::make_pair(1.0, 0.0); // wind direction

	std::vector<std::pair<double, double>> initial_h; // initial height fields for each K
    std::vector<std::pair<double, double>> current_h; // current height fields for each K


    const double D = 6.0; // Depth below mean water level (for dispersion relation)


    std::pair<double, double> index_to_k_vector
	(
		int k_index
	);

    std::pair<double, double> index_to_horizontal_pos(int k_index);

	std::pair<double, double> sample_complex_gaussian
		();
    double phillips_spectrum(std::pair<double, double> k);
    std::pair<double, double> amplitude_0(std::pair<double, double> k);
	int k_index_to_negative_k_index
	(
		int k_index
	);
	double omega_dispersion
	(
		double k_magnitude,
		bool is_shallow=false
	);
    std::pair<double, double> exp_complex(double ix);

    std::pair<double, double> amplitude_t(double t, int index, std::pair<double, double> k);

};


#endif //OCEAN_H
