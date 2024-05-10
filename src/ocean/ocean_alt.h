#ifndef OCEAN_ALT_H
#define OCEAN_ALT_H
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#define EIGEN_DONT_VECTORIZE

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

struct FoamConstants{
    std::vector<Eigen::Vector3f> positions;
    std::vector<Eigen::Vector2f> k_vectors;
    std::vector<float> wavelengths;
    std::vector<Eigen::Vector2f> texCoords;
};

struct OceanSpray{
    Eigen::Vector3f height;
    Eigen::Vector3f slope;

};

class ocean_alt
{
public:
    ocean_alt();
    void updateVertexAmplitudes(double t);
	void update_ocean();
    std::vector<Eigen::Vector3f> get_vertices();
    std::vector<Eigen::Vector3i> get_faces();
    void fft_prime(double t);
    std::vector<Eigen::Vector3f> getNormals();

    FoamConstants getFoamConstants(){
        return m_foam_constants;
    }

    std::vector<Eigen::Vector3f> m_vertices; // current displacement vector for each K
    std::vector<OceanSpray> m_heights; // stores height above threshold








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

    // FOAM
    std::vector<float> m_saturation;

    std::map<int, WaveIndexConstant> m_waveIndexConstants; // stores constants that only need to be calculate once for each grid constant

	const double Lx = 512.0;
	const double Lz = 512.0;

    const int num_rows = 128;
    const int num_cols = 128;

    const int num_tiles_x = 1;
    const int num_tiles_z = 1;

	const double vertex_displacement = Lx / 2;

	const int N = num_rows*num_cols; // total number of grid points
	const double lambda = .5; // how much displacement matters
	const double spacing = 1.0; // spacing between grid points

    const double A = 100; // numeric constant for the Phillips spectrum
    const double V = 100; // wind speed
    const double gravity = 9.81;
    const double L = V*V/gravity;
    const Eigen::Vector2d omega_wind = Eigen::Vector2d(1.0, 0.0); // wind direction, used in Phillips equation

    std::vector<Eigen::Vector2d> m_current_h; // current height fields for each K
   	// std::vector<Eigen::Vector2d> m_displacements; // current displacement vector for each K
    // std::vector<Eigen::Vector2d> m_slopes; // current displacement vector for each K
	std::vector<Eigen::Vector2d> m_slopes_x;
	std::vector<Eigen::Vector2d> m_slopes_z;
	std::vector<Eigen::Vector2d> m_displacements_x;
	std::vector<Eigen::Vector2d> m_displacements_z;
    //std::vector<Eigen::Vector3f> m_slope_vectors; // current displacement vector for each K

    std::vector<Eigen::Vector3f> m_normals; // normal calculations

    // FOR FOAM:
    FoamConstants m_foam_constants;
    float height_threshold = 2.f;



    float max = 0;
    float min = 0;
    int iterations = 0;



    const double D = 1.0; // Depth below mean water level (for dispersion relation)


	std::vector<Eigen::Vector2d> fast_fft
	(
		std::vector<Eigen::Vector2d> h
	);
};

#endif // OCEAN_ALT_H
