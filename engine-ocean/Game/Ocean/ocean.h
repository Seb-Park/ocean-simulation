#ifndef ocean_H
#define ocean_H

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

class ocean
{
public:
    ocean();
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
    std::pair<double, double> sample_complex_gaussian();
    std::vector<Eigen::Vector2d> fast_fft(std::vector<Eigen::Vector2d> h);









    std::map<int, WaveIndexConstant> m_waveIndexConstants; // stores constants that only need to be calculate once for each grid constant



    const double Lx = 10.0;
    const double Lz = 10.0;

    const int num_rows = 8;
    const int num_cols = 8;

    const int N = num_rows*num_cols; // total number of grid points
    const double lambda = .40; // how much displacement matters
    const double spacing = 35.0; // spacing between grid points

    const double A = 1.0; // numeric constant for the Phillips spectrum
    const double V = 5.5; // wind speed
    const double gravity = 9.81;
    const double L = V*V/gravity;
    const Eigen::Vector2d omega_wind = Eigen::Vector2d(1.0, 0.0); // wind direction, used in Phillips equation

    std::vector<Eigen::Vector2d> m_current_h; // current height fields for each K
    std::vector<Eigen::Vector2d> m_displacements; // current displacement vector for each K



    const double D = 1.0; // Depth below mean water level (for dispersion relation)


};

#endif // ocean_H