#include "ocean_alt.h"
#include <iostream>


ocean_alt::ocean_alt()
{
    // to be used for efficiency during fft
    std::cout << "hello" << std::endl;
    init_wave_index_constants();

}

// initializes static constants (aka they are not time dependent)
void ocean_alt::init_wave_index_constants(){

    for (int i=0; i<N; i++){
        Eigen::Vector2i m_n = index_1d_to_2d(i);
        int n_prime = m_n[0];
        int m_prime = m_n[1];

        Eigen::Vector2d k = get_k_vector(n_prime, m_prime);
        Eigen::Vector2d k_conj = get_k_vector(-n_prime, m_prime);


        // store h0'(n,m) and w'(n,m) for every index, to be used for later
        Eigen::Vector2d h0_prime = h_0_prime(k);

        // conjugate of a+bi is a-bi
        Eigen::Vector2d h0_prime_conj = h_0_prime(k_conj);
        h0_prime_conj = Eigen::Vector2d(h0_prime_conj[0], -h0_prime_conj[1]);

        double w_prime = omega_prime(k);

        // populate map to be used for later
        WaveIndexConstant wave_const;
        wave_const.h0_prime = h0_prime;
        wave_const.h0_prime_conj = h0_prime_conj;
        wave_const.w_prime = w_prime;
        wave_const.w_prime = w_prime;
        wave_const.base_horiz_pos = get_horiz_pos(i);
        wave_const.k_vector = k;

        m_waveIndexConstants[i] = wave_const;

        // initialize m_current_h to be h0 for now
        m_current_h.push_back(h0_prime);
        m_displacements.push_back(Eigen::Vector2d(0.0, 0.0));
        m_slopes.push_back(Eigen::Vector2d(0.0, 0.0));
        m_normals.push_back(Eigen::Vector3f(0.0, 1.0, 0.0));

    }
}

// fast fourier transform at time t
void ocean_alt::fft_prime(double t){

    // FFT
    std::vector<Eigen::Vector2d> h_tildas = std::vector<Eigen::Vector2d>();

    // find each h_tilda at each index, to be used for next for loop
    for (int i=0; i<N; i++){
        Eigen::Vector2d h_t_prime = h_prime_t(i, t); // vector(real, imag)

        h_tildas.emplace_back(h_t_prime);
    }

	bool fast = false;
	if (fast)
	{
		std::vector<Eigen::Vector2d> tmp = fast_fft(h_tildas);
		for (int i = 0; i < N; i++)
		{
			m_current_h[i] = tmp[i];

			// update displacements and slopes
			Eigen::Vector2d k_vector = m_waveIndexConstants[i].k_vector;
			Eigen::Vector2d k_normalized = k_vector.normalized();

			double imag_comp = m_current_h[i][1];
			m_displacements[i] += k_normalized*imag_comp;
			m_slopes[i] += k_vector*imag_comp;
		}

		return;
	}

    // for each position in grid, sum up amplitudes dependng on that position
    for (int i=0; i<N; i++){
        Eigen::Vector2d x_vector = m_waveIndexConstants[i].base_horiz_pos;
        m_current_h[i] = Eigen::Vector2d(0.0, 0.0);
        m_displacements[i] = Eigen::Vector2d(0.0, 0.0);
        m_slopes[i] = Eigen::Vector2d(0.0, 0.0);



        for (int j = 0; j < N; j++){
            Eigen::Vector2d k_vector = m_waveIndexConstants[j].k_vector;
            Eigen::Vector2d h_tilda_prime = h_tildas[j]; // vector(real, imag)


            // add x vector and k vector as imaginary numbers
            double imag_xk_sum = x_vector.dot(k_vector);
            Eigen::Vector2d exp = complex_exp(-imag_xk_sum); // vector(real, imag)

            double real_comp = h_tilda_prime[0]*exp[0] - h_tilda_prime[1]*exp[1];
            double imag_comp = h_tilda_prime[0]*exp[1] + h_tilda_prime[1]*exp[0];

            m_current_h[i] += Eigen::Vector2d(real_comp, imag_comp);

            Eigen::Vector2d k_normalized = k_vector.normalized();

            m_displacements[i] += k_normalized*imag_comp;
            m_slopes[i] += -k_vector*imag_comp;
        }
    }

}

// time dependent calculation of h'(n,m,t)
Eigen::Vector2d ocean_alt::h_prime_t(int i, double t){
    Eigen::Vector2d h0_prime = m_waveIndexConstants[i].h0_prime; // vector(real, imag)
    Eigen::Vector2d h0_prime_conj = m_waveIndexConstants[i].h0_prime_conj; // vector(real, imag)
    double w_prime = m_waveIndexConstants[i].w_prime;

    Eigen::Vector2d pos_complex_exp = complex_exp(w_prime*t); // vector(real, imag)
    Eigen::Vector2d neg_complex_exp = complex_exp(-w_prime*t); // vector(real, imag)

    // now multiply our four vector(real, imag) out

    double real_comp =
            h0_prime[0]*pos_complex_exp[0]
            - h0_prime[1]*pos_complex_exp[1]
            + h0_prime_conj[0]*neg_complex_exp[0]
            + h0_prime_conj[1]*neg_complex_exp[1];

    double imag_comp =
            h0_prime[0]*pos_complex_exp[1]
            + h0_prime[1]*pos_complex_exp[0]
            + h0_prime_conj[0]*neg_complex_exp[1]
            - h0_prime_conj[1]*neg_complex_exp[0];



    return Eigen::Vector2d(real_comp, imag_comp);
}

double ocean_alt::omega_prime(Eigen::Vector2d k){
    // calculate omega^4 first to prevent sqrts
    double w = sqrt(gravity*k.norm());

    return w;
}

Eigen::Vector2d ocean_alt::h_0_prime(Eigen::Vector2d k){
    double Ph_prime = phillips_prime(k);
    std::pair<double,double> randoms = sample_complex_gaussian();
    double random_r = randoms.first;
    double random_i = randoms.second;

    // seperate real and imag products
    double coeff = 0.707106781187 * sqrt(Ph_prime);
    double real_comp = coeff*random_r;
    double imag_comp = coeff*random_i;

    return Eigen::Vector2d(real_comp, imag_comp);
}


double ocean_alt::phillips_prime(Eigen::Vector2d k){
    double k_mag = k.norm();

    k.normalize();
    double dot_prod = k.dot(omega_wind);

    double output = 0.0;
    // l = 1
    if (k_mag < .0001) return 0.0;

    if (k_mag > 1.0){

       output =  A*exp(-(k_mag*k_mag))*dot_prod*dot_prod/(k_mag*k_mag*k_mag*k_mag);
    } else {
       output =  A*exp(-1.0/(k_mag*L*k_mag*L))*dot_prod*dot_prod/(k_mag*k_mag*k_mag*k_mag);

    }



    return output;
}

Eigen::Vector2d ocean_alt::get_k_vector(int n_prime, int m_prime){
    double n_ = (double)n_prime;
    double m_ = (double)m_prime;
    double N_ = (double)num_rows;
    double M_ = (double)num_cols;

    double k_x = (2.0*M_PI*n_ - M_PI*N_)/Lx;
    double k_z = (2.0*M_PI*m_ - M_PI*M_)/Lz;

    return Eigen::Vector2d(k_x, k_z);
}

Eigen::Vector2d ocean_alt::get_horiz_pos(int i){
    Eigen::Vector2i m_n = index_1d_to_2d(i);
    double n_prime = (double)m_n[0];
    double m_prime = (double)m_n[1];
    double N_ = (double)num_rows;
    double M_ = (double)num_cols;


    double x = (n_prime-.5*N_)*Lx / N_;
    double z = (m_prime-.5*M_)*Lz / M_;



    return Eigen::Vector2d(x, z);
}


Eigen::Vector2i ocean_alt::index_1d_to_2d(int i){
    int row = i/num_rows; // n'
    int col = i%num_rows; // m'

    return Eigen::Vector2i(row, col);

}

std::pair<double,double> ocean_alt::sample_complex_gaussian(){
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

Eigen::Vector2d ocean_alt::complex_exp(double exponent){
    double real = cos(exponent);
    double imag = sin(exponent);

    return Eigen::Vector2d(real, imag);
}

std::vector<Eigen::Vector3f> ocean_alt::get_vertices()
{
    std::vector<Eigen::Vector3f> vertices = std::vector<Eigen::Vector3f>();
    for (int i = 0; i < N; i++){
        Eigen::Vector2d horiz_pos = spacing*m_waveIndexConstants[i].base_horiz_pos;
        Eigen::Vector2d amplitude = m_current_h[i];
        float height = amplitude[0];

        Eigen::Vector2d slope = m_slopes[i] * .3f;
        Eigen::Vector3f s = Eigen::Vector3f(-slope[0], 0.0, -slope[1]);
        Eigen::Vector3f y = Eigen::Vector3f(0.0, 1.0, 0.0);

        float xs = 1.f + s[0]*s[0];
        float ys = 1.f + s[1]*s[1];
        float zs = 1.f + s[2]*s[2];

        Eigen::Vector3f diff = y - s;
        // Eigen::Vector3f norm = Eigen::Vector3f(diff[0]/ sqrt(xs), diff[1]/ sqrt(ys), diff[2]/sqrt(zs));

        // NEW
        Eigen::Vector3f norm = Eigen::Vector3f(-slope[0], 1.0, -slope[1]);
        norm.normalize();
        // NEW





        //if (i==6) std::cout << amplitude[0] << std::endl;

        // calculate displacement
        Eigen::Vector2d disp = lambda*m_displacements[i];

        //


        // for final vertex position, use the real number component of amplitude vector
        vertices.push_back(Eigen::Vector3f(horiz_pos[0] + disp[0], height, horiz_pos[1] + disp[1]));
        m_normals[i] = norm.normalized();//Eigen::Vector3f(-slope[0], 1.0, -slope[1]).normalized();
        //std::cout << "normal: " << m_normals[i] << std::endl;

    }
    return vertices;
}

std::vector<Eigen::Vector3f> ocean_alt::getNormals(){
    return m_normals;
}

std::vector<Eigen::Vector3i> ocean_alt::get_faces()
{
    // connect the vertices into faces
    std::vector<Eigen::Vector3i> faces = std::vector<Eigen::Vector3i>();
    for (int i = 0; i < N; i++)
    {
        int x = i / num_rows;
        int z = i % num_rows;

        // connect the vertices into faces
        if (x < num_rows - 1 && z < num_cols - 1)
        {
            int i1 = i;
            int i2 = i + 1;
            int i3 = i + num_rows;
            int i4 = i + num_rows + 1;

            faces.emplace_back(i2, i1, i3);
            faces.emplace_back(i2, i3, i4);
                        faces.emplace_back(i1, i2, i3);
                        faces.emplace_back(i3, i2, i4);
        }
    }
    return faces;
}

std::vector<Eigen::Vector2d> ocean_alt::fast_fft
	(
		std::vector<Eigen::Vector2d> h
	)
{
	int N = h.size();
	int exponent = 0;
	int power_of_2 = 1;
	while (power_of_2 < N)
	{
		power_of_2 *= 2;
		exponent++;
	}

	std::vector<Eigen::Vector2d> H = std::vector<Eigen::Vector2d>(N);
	// pad with zeros
	for (int i = 0; i < N; i++)
	{
		H[i] = h[i];
	}
	for (int i = N; i < power_of_2; i++)
	{
		H.emplace_back(Eigen::Vector2d(0.0, 0.0));
	}

	// bit reverse the indices of the input data array
	std::vector<Eigen::Vector2d> temp = std::vector<Eigen::Vector2d>(power_of_2);
	for (int i = 0; i < power_of_2; i++)
	{
		int j = 0;
		for (int k = 0; k < exponent; k++)
		{
			j = (j << 1) | (i >> k & 1);
		}
		temp[j] = H[i];
	}
	for (int i = 0; i < power_of_2; i++)
	{
		H[i] = temp[i];
	}

	// outer loop through levels i->p of even-odd decompositions, beginning with the final decompositions
	// where individual y_m are used through the first decomposition where Y_n^e and Y_n^o are involved
	for (int i = 1; i <= exponent; i++)
	{

		// nested middle loop where each level is grouped into terms according to values of k in the factor
		// W_N^k = e^(-2*pi*i*k/N) appearing in each sum
		int N_over_2 = 1 << i;
		int N_over_4 = 1 << (i - 1);
		double theta = 2 * M_PI / N_over_2;
		Eigen::Vector2d W_N = Eigen::Vector2d(cos(theta), -sin(theta));
		Eigen::Vector2d W = Eigen::Vector2d(1.0, 0.0);
		for (int k = 0; k < N_over_4; k++)
		{
			// innermost loop where each group is split into individual sums
			for (int j = k; j < power_of_2; j += N_over_2)
			{
				Eigen::Vector2d U = H[j];
				Eigen::Vector2d V = H[j + N_over_4];
				H[j] = U + W[0] * V - W[1] * V;
				H[j + N_over_4] = U - (W[0] * V - W[1] * V);
			}
			W = W[0] * W_N - W[1] * W_N;
		}
	}

	return H;
}
