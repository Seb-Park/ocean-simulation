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
    float tex_step = 1.f/num_rows;

    for (int i=0; i<N; i++){
        Eigen::Vector2i m_n = index_1d_to_2d(i);
        int n_prime = m_n[0];
        int m_prime = m_n[1];

        Eigen::Vector2d k = get_k_vector(n_prime, m_prime);
        Eigen::Vector2d k_conj = get_k_vector(-n_prime, m_prime);


//        Eigen::Vector3f v = Eigen::Vector3f(0,0,1);
//        Eigen::Vector3f norm = Eigen::Vector3f(0,1,0);
//        if (abs(norm[1]) < 1.f){
//            v = (Eigen::Vector3f(0,1,0) - norm[1]*norm).normalized();
//        }
//        Eigen::Vector3f u = norm.cross(v).normalized();

//        float u_coord = u.dot(Eigen::Vector3f(n_prime, 0, m_prime)) / 64.f;
//        float v_coord = v.dot(Eigen::Vector3f(n_prime, 0, m_prime)) / 64.f;

//        //std::cout << u_coord << ", " << v_coord << std::endl;


        // texture coord:
        Eigen::Vector2f texCoord = Eigen::Vector2f(1, 1);


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
        wave_const.base_horiz_pos = get_horiz_pos(i);
        wave_const.k_vector = k;

        m_waveIndexConstants[i] = wave_const;

        // initialize m_current_h to be h0 for now
//        m_current_h.push_back(h0_prime);
        m_current_h.push_back(Eigen::Vector2d(0.0, 0.0));
        // m_displacements.push_back(Eigen::Vector2d(0.0, 0.0));
        // m_slopes.push_back(Eigen::Vector2d(0.0, 0.0));
        m_normals.push_back(Eigen::Vector3f(0.0, 1.0, 0.0));

        // initialize foam constant vectors
        m_foam_constants.k_vectors.push_back(Eigen::Vector2f(k[0], k[1]));
        m_foam_constants.positions.push_back(Eigen::Vector3f(0,0,0));
        m_foam_constants.wavelengths.push_back(0);
       // m_foam_constants.texCoords.push_back(texCoord);


		m_slopes_x.push_back(Eigen::Vector2d(0.0, 0.0));
		m_slopes_z.push_back(Eigen::Vector2d(0.0, 0.0));

		m_displacements_x.push_back(Eigen::Vector2d(0.0, 0.0));
		m_displacements_z.push_back(Eigen::Vector2d(0.0, 0.0));

		m_vertices.push_back(Eigen::Vector3f(0.0, 0.0, 0.0));
    }
}

// fast fourier transform at time t
void ocean_alt::fft_prime(double t){

    // FFT
	std::vector<Eigen::Vector2d> h_tildas = std::vector<Eigen::Vector2d>();
	// std::vector<Eigen::Vector2d> ikh = std::vector<Eigen::Vector2d>();
//	std::vector<Eigen::Vector2d> neg_ik_hat_h = std::vector<Eigen::Vector2d>();

	std::vector<Eigen::Vector2d> ikhx = std::vector<Eigen::Vector2d>();
	std::vector<Eigen::Vector2d> ikhz = std::vector<Eigen::Vector2d>();

	std::vector<Eigen::Vector2d> neg_ik_hat_h_x = std::vector<Eigen::Vector2d>();
	std::vector<Eigen::Vector2d> neg_ik_hat_h_z = std::vector<Eigen::Vector2d>();

    // find each h_tilda at each index, to be used for next for loop
	for (int i=0; i<N; i++){
		Eigen::Vector2d h_t_prime = h_prime_t(i, t); // vector(real, imag)

		h_tildas.emplace_back(h_t_prime);

		Eigen::Vector2d k_vector = m_waveIndexConstants[i].k_vector;
		// ikh.emplace_back(-h_t_prime[1] * k_vector[0], -h_t_prime[1] * k_vector[1]);
		ikhx.emplace_back(-h_t_prime[1] * k_vector[0], -h_t_prime[1] * k_vector[0]);
		ikhz.emplace_back(-h_t_prime[1] * k_vector[1], -h_t_prime[1] * k_vector[1]);


//		Eigen::Vector2d neg_ik_hat_h_val =
//			Eigen::Vector2d(k_normalized[1] * h_t_prime[1], k_normalized[0] * h_t_prime[1]);
//		neg_ik_hat_h.emplace_back(neg_ik_hat_h_val);
		double len = k_vector.norm();
		if (len < .000001)
		{
			neg_ik_hat_h_x.emplace_back(0.0, 0.0);
			neg_ik_hat_h_z.emplace_back(0.0, 0.0);
		}
		else
		{
			Eigen::Vector2d k_normalized = k_vector.normalized();
			neg_ik_hat_h_x.emplace_back(k_normalized[0] * h_t_prime[1], k_normalized[0] * h_t_prime[1]);
			neg_ik_hat_h_z.emplace_back(k_normalized[1] * h_t_prime[1], k_normalized[1] * h_t_prime[1]);
		}
	}

	bool fast = true;
	if (fast)
	{
		std::vector<Eigen::Vector2d> tmp = fast_fft(h_tildas);
//		std::vector<Eigen::Vector2d> tmp2 = fast_fft(ikh);
		// std::vector<Eigen::Vector2d> tmp3 = fast_fft(neg_ik_hat_h);
		std::vector<Eigen::Vector2d> tmp4 = fast_fft(ikhx);
		std::vector<Eigen::Vector2d> tmp5 = fast_fft(ikhz);
		std::vector<Eigen::Vector2d> tmp6 = fast_fft(neg_ik_hat_h_x);
		std::vector<Eigen::Vector2d> tmp7 = fast_fft(neg_ik_hat_h_z);
		for (int i = 0; i < N; i++)
		{
			m_current_h[i] = tmp[i];
			// m_slopes[i] = tmp2[i];
			// m_displacements[i] = tmp3[i];
			m_slopes_x[i] = tmp4[i];
			m_slopes_z[i] = tmp5[i];
			m_displacements_x[i] = tmp6[i];
			m_displacements_z[i] = tmp7[i];
		}

		return;
	}

    // for each position in grid, sum up amplitudes dependng on that position
    for (int i=0; i<N; i++){
        Eigen::Vector2d x_vector = m_waveIndexConstants[i].base_horiz_pos;
        m_current_h[i] = Eigen::Vector2d(0.0, 0.0);
        // m_displacements[i] = Eigen::Vector2d(0.0, 0.0);
        // m_slopes[i] = Eigen::Vector2d(0.0, 0.0);



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

            // m_displacements[i] += k_normalized*imag_comp;
            //  m_slopes[i] += -k_vector*imag_comp;
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
		(h0_prime[0]*pos_complex_exp[0] - h0_prime[1]*pos_complex_exp[1]) +
		(h0_prime_conj[0]*neg_complex_exp[0] - h0_prime_conj[1]*neg_complex_exp[1]);

    double imag_comp =
		(h0_prime[0]*pos_complex_exp[1] + h0_prime[1]*pos_complex_exp[0]) +
		(h0_prime_conj[0]*neg_complex_exp[1] + h0_prime_conj[1]*neg_complex_exp[0]);

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

    double k_x = (2*M_PI*n_ - M_PI*N_)/Lx;
    double k_z = (2*M_PI*m_ - M_PI*M_)/Lz;

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
	double real = sqrt(-2*log(uniform_1)) * cos(2*M_PI*uniform_2);
	double imag = sqrt(-2*log(uniform_1)) * sin(2*M_PI*uniform_2);

    return std::make_pair(real, imag);
}

Eigen::Vector2d ocean_alt::complex_exp(double exponent){
    double real = cos(exponent);
    double imag = sin(exponent);

    return Eigen::Vector2d(real, imag);
}

void ocean_alt::update_ocean()
{
    std::vector<Eigen::Vector3f> vertices = std::vector<Eigen::Vector3f>();
	// reset normals & vertices arrays for the single tile
    m_vertices = std::vector<Eigen::Vector3f>(N);
	m_normals = std::vector<Eigen::Vector3f>(N);
    for (int i = 0; i < N; i++){
        Eigen::Vector2d horiz_pos = spacing*m_waveIndexConstants[i].base_horiz_pos;
        Eigen::Vector2d amplitude = m_current_h[i];
        float height = amplitude[0];

		if (iterations++ > 1)
		{
			if (height < min) min = height;
			if (height > max)
			{
				max = height;
//				std::cout << "changed!! max: " << max << std::endl;
			}
		}


		// Eigen::Vector2d slope = m_slopes[i] * .3f;
        // Eigen::Vector3f s = Eigen::Vector3f(-slope[0], 0.0, -slope[1]);
        // Eigen::Vector3f y = Eigen::Vector3f(0.0, 1.0, 0.0);

//        float xs = 1.f + s[0]*s[0];
//        float ys = 1.f + s[1]*s[1];
//        float zs = 1.f + s[2]*s[2];
//
//        Eigen::Vector3f diff = y - s;
//        Eigen::Vector3f Eigen::Vector3f(diff[0]/ sqrt(xs), diff[1]/ sqrt(ys), diff[2]/sqrt(zs));

         // NEW
		 Eigen::Vector3f norm = Eigen::Vector3f(-m_slopes_x[i][0], 1.0, -m_slopes_z[i][0]);
         norm = norm.normalized(); // FIXME: why do I have to be inverted?

        //if (i==6) std::cout << amplitude[0] << std::endl;

        // calculate displacement
        // Eigen::Vector2d disp = lambda*m_displacements[i];
		Eigen::Vector2d disp = lambda*Eigen::Vector2d(m_displacements_x[i][0], m_displacements_z[i][0])
			+ Eigen::Vector2d(vertex_displacement, vertex_displacement); // set corner at 0,0 for retiling


        // for final vertex position, use the real number component of amplitude vector
		m_vertices[i] = {horiz_pos[0] + disp[0], height, horiz_pos[1] + disp[1]};
        m_normals[i] = norm.normalized();//Eigen::Vector3f(-slope[0], 1.0, -slope[1]).normalized();
        //std::cout << "normal: " << m_normals[i] << std::endl
        Eigen::Vector2i m_n = index_1d_to_2d(i);


       // m_foam_constants.wavelengths[i] = 2.f* M_PI * m_slopes[i].dot(m_slopes[i]) /  Lx;
        float h_0 = 0; // min*.2f;
        float h_max = max*.35f; // the smaller the constant, the more foam there is
        m_foam_constants.wavelengths[i] = (height - h_0 ) / (h_max - h_0);

//        if (i < 5){
//            std::cout << h_0 << ", " << h_max << std::endl;
//            std::cout << m_foam_constants.wavelengths[i] << std::endl;
//        }



    }


    // populate foam constants
    m_foam_constants.positions = vertices;
}

std::vector<Eigen::Vector3f> ocean_alt::get_vertices(){
	// extend the returned array based on the tilecount
	std::vector<Eigen::Vector3f> vertices = std::vector<Eigen::Vector3f>();
	for (int i = 0; i < num_tiles_x; i++)
	{
		for (int j = 0; j < num_tiles_z; j++)
		{
			for (int k = 0; k < N; k++)
			{
				double c = Lx - 2 / (num_rows / Lx);
				Eigen::Vector3f vertex = m_vertices[k] + Eigen::Vector3f(i*(c), 0.0, (j)*(c));
				vertices.push_back(vertex);
			}
		}
	}

	return vertices;
}

std::vector<Eigen::Vector3f> ocean_alt::getNormals(){
	// based on the tile count, add more to the normals
	std::vector<Eigen::Vector3f> normals = std::vector<Eigen::Vector3f>();
	// do the x 1D direction first
	for (int i = 0; i < num_tiles_x; i++)
	{
		for (int j = 0; j < num_tiles_z; j++)
		{
			for (int k = 0; k < N; k++)
			{
				normals.push_back(m_normals[k]);
			}
		}
	}


    return normals;
}

std::vector<Eigen::Vector3i> ocean_alt::get_faces()
{
    // connect the vertices into faces
    std::vector<Eigen::Vector3i> faces = std::vector<Eigen::Vector3i>();
//	for (int i = 0; i < num_tiles_x; i++)
//	{
//		for (int j = 0; j < num_tiles_z; j++)
//		{
//			for (int k = 0; k < N; k++)
//			{
//				int x = k % num_rows;
//				int z = k / num_rows;
//
//				// connect the vertices into faces
//				if (x < num_rows - 1 && z < num_cols - 1)
//				{
//					int tile_index_offset = (j + num_tiles_z * i) * N;
//					int i1 = k + tile_index_offset;
//					int i2 = k + 1 + tile_index_offset;
//					int i3 = k + num_rows + tile_index_offset;
//					int i4 = k + num_rows + 1 + tile_index_offset;
//
//					faces.emplace_back(i2, i1, i3);
//					faces.emplace_back(i2, i3, i4);
//				}
//			}
//		}
//	}
//
//	return faces;


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

Eigen::Vector2d muliply_complex(Eigen::Vector2d a, Eigen::Vector2d b)
{
	double real = a[0] * b[0] - a[1] * b[1];
	double imag = a[0] * b[1] + a[1] * b[0];
	return Eigen::Vector2d(real, imag);
}

std::vector<Eigen::Vector2d> ifft_1d
	(
		std::vector<Eigen::Vector2d> frequencies,
		bool is_vertical
	)
{
	// one D case, assuming is a square
	int N = frequencies.size();
	// make two buffers for intermediate butterfly values
	std::vector<Eigen::Vector2d> buffer1 = std::vector<Eigen::Vector2d>(N);
	std::vector<Eigen::Vector2d> buffer2 = std::vector<Eigen::Vector2d>(N);

	// fill buffer one with the frequencies in bit reverse order
	int log2_N = log2(N);
	for (int i = 0; i < N; i++)
	{
		int reversed = 0;
		for (int j = 0; j < log2_N; j++)
		{
			reversed |= ((i >> j) & 1) << (log2_N - 1 - j);
		}
		// std::cout << "reversed, i: " << reversed << ", " << i << std::endl;
		buffer1[i] = frequencies[reversed];
	}
	bool reading_buffer1 = true;

	// go over the stages
	for (int stage = 1; stage <= log2_N; stage++)
	{
		// go over the groups
		for (int group = 0; group < N / pow(2, stage); group++)
		{
			// go over the butterflies
			for (int butterfly = 0; butterfly < pow(2, stage - 1); butterfly++)
			{
				// calculate the indices
				int index1 = group * pow(2, stage) + butterfly;
				int index2 = group * pow(2, stage) + pow(2, stage - 1) + butterfly;

				// calculate the twiddle factor
				int index = group * pow(2, stage) + butterfly;
				float w = index * 2 * M_PI / pow(2, stage);
				Eigen::Vector2d twiddle_factor = {cos(w), sin(w)};

				if (reading_buffer1)
				{
					buffer2[index1] = buffer1[index1] + muliply_complex(twiddle_factor, buffer1[index2]);
					buffer2[index2] = buffer1[index1] - muliply_complex(twiddle_factor, buffer1[index2]);
				}
				else
				{
					buffer1[index1] = buffer2[index1] + muliply_complex(twiddle_factor, buffer2[index2]);
					buffer1[index2] = buffer2[index1] - muliply_complex(twiddle_factor, buffer2[index2]);
				}
			}
		}
		reading_buffer1 = !reading_buffer1;
	}

	// return the buffer that was read last
	if (reading_buffer1)
	{
		return buffer1;
	}
	else
	{
		return buffer2;
	}
}


std::vector<Eigen::Vector2d> ocean_alt::fast_fft
	(
		std::vector<Eigen::Vector2d> h
	)
{
	// do a vertical fft on each column
	for (int i = 0; i < num_rows; i++)
	{
		std::vector<Eigen::Vector2d> col = std::vector<Eigen::Vector2d>();
		for (int j = 0; j < num_cols; j++)
		{
			col.push_back(h[i + j * num_rows]);
		}
		std::vector<Eigen::Vector2d> col_fft = ifft_1d(col, true);
		for (int j = 0; j < num_cols; j++)
		{
			h[i + j * num_rows] = col_fft[j];
		}
	}

	// do a horizontal fft on each row
	for (int i = 0; i < num_cols; i++)
	{
		std::vector<Eigen::Vector2d> row = std::vector<Eigen::Vector2d>();
		for (int j = 0; j < num_rows; j++)
		{
			row.push_back(h[i * num_rows + j]);
		}
		std::vector<Eigen::Vector2d> row_fft = ifft_1d(row, false);
		for (int j = 0; j < num_rows; j++)
		{
			h[i * num_rows + j] = row_fft[j];
		}
	}

	// divide by N*N and add the signs based on the indices
    double sign[] = {-1.0, 1.0};
	for (int i = 0; i < N; i++)
	{
		h[i] /= N;
		// h[i] /= sqrt(N);
		h[i] *= sign[(i / num_rows + i % num_cols) % 2];
	}

	return h;
}
