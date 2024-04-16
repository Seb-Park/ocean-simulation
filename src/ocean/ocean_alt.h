#ifndef OCEAN_ALT_H
#define OCEAN_ALT_H

#include <map>
#include <vector>
#include <utility>
#include <Eigen/Dense>

// for every 1d index up to length*width
struct WaveIndexConstant{
    Eigen::Vector2d h0_prime = Eigen::Vector2d(0.f, 0.f);
    Eigen::Vector2d h0_prime_conj = Eigen::Vector2d(0.f, 0.f);

    double w_prime = 0.0;


    Eigen::Vector2d base_horiz_pos = Eigen::Vector2d(0.f, 0.f); // static horiz pos with no displacement
    Eigen::Vector2d k_vector = Eigen::Vector2d(0.f, 0.f); // static horiz pos with no displacement



};

class ocean_alt
{
public:
    ocean_alt();
    void updateVertexAmplitudes(double t);
    std::vector<Eigen::Vector3f> get_vertices();
    std::vector<Eigen::Vector3i> get_faces();
    void fft_prime(double t);





private:

    Eigen::Vector2i index_1d_to_2d(int i);
    Eigen::Vector2d get_k_vector(int n_prime, int m_prime);
    double phillips_prime(Eigen::Vector2d k);
    Eigen::Vector2d h_0_prime(Eigen::Vector2d k);
    double omega_prime(Eigen::Vector2d k);
    void init_wave_index_constants();
    Eigen::Vector2d complex_exp(double exponent);
    Eigen::Vector2d h_prime_t(int i, double t);
    Eigen::Vector2d get_horiz_pos(int i);







    std::map<int, WaveIndexConstant> m_waveIndexConstants; // stores constants that only need to be calculate once for each grid constant




    const int length = 10; // length of grid (L_x)
    const int width = 10; // width of grid (L_z)
    const int N = length * width; // total number of grid points
    const double spacing = 1.0; // spacing between grid points
    const double lambda = 4.0; // how much displacement matters

    const double A = 1.0; // numeric constant for the Phillips spectrum
    const double V = 1.0; // wind speed
    const double gravity = 9.8;
    const double L = V*V/gravity;
    const Eigen::Vector2d omega_wind = Eigen::Vector2d(1.0, 0.0); // wind direction, used in Phillips equation

    std::vector<std::pair<double, double>> initial_h; // initial height fields for each K
    std::vector<Eigen::Vector2d> m_current_h; // current height fields for each K


    const double D = 1.0; // Depth below mean water level (for dispersion relation)


    std::pair<double, double> k_index_to_k_vector
    (
        int k_index
    );

    std::pair<double, double> k_index_to_horiz_pos(int k_index);

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

    std::pair<double, double> calculate_displacement(int k_index);

};

#endif // OCEAN_ALT_H
