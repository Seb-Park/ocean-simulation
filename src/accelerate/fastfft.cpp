//
// Created by Michael Foiani on 4/29/24.
//

#include "fastfft.h"
#include <iostream>

// constructor
fastfft::fastfft(int N)
{
	// assert N is a power of 2
	assert((N & (N - 1)) == 0);

	// set variables assisvaited to N
	this->N = N;
	this->log_2_N = (unsigned int) log2(N);

	// initialize reverse bit indices
	reversed_bit_indices = std::vector<int>(N);
	for (int i = 0; i < N; i++)
	{
		int tmp = i;
		int res = 0;
		for (int j = 0; j < log_2_N; j++) {
			res = (res << 1) + (tmp & 1);
			tmp >>= 1;
		}
		reversed_bit_indices[i] = res;
	}

	// fill in the butterfly map
	for (int stage = 0; stage < log_2_N; stage++)
	{
		m_butterfly_map.push_back(std::vector<Eigen::Vector4d>(N));
		for (int n = 0; n < N; n++)
		{
			precompute_butterfly_map(stage, n);
		}
	}


	int pow2 = 1;
	m_butterfly = std::vector<std::vector<complex>>(log_2_N);
	for (int stage = 0; stage < log_2_N; stage++)
	{
		m_butterfly[stage] = std::vector<complex>(pow2);
		for (int n = 0; n < pow2; n++)
		{
			m_butterfly[stage][n] =
				{
				cos(M_PI * n / pow2), sin(M_PI * n / pow2)
				};
		};
		pow2 *= 2;
	}
	c = std::vector<std::vector<complex>>(2);
	c[0] = std::vector<complex>(N);
	c[1] = std::vector<complex>(N);
	which = 0;

	// inialize variables for the computation
	m_pingpong = 0;
	m_pingpong0 = std::vector<std::vector<Eigen::Vector2d>>(N, std::vector<Eigen::Vector2d>(N));
	m_pingpong1 = std::vector<std::vector<Eigen::Vector2d>>(N, std::vector<Eigen::Vector2d>(N));

	// initialize results
	results = std::vector<std::vector<complex>>(N, std::vector<complex>(N));
}

void fastfft::precompute_butterfly_map
(
	int stage,
	int n
)
{
	// calculate the twiddle from k
	double k = ((int) (((double) N) / pow(2, stage + 1)) * n) % N;
	double t = 2 * M_PI * k / (double) N;
	complex twiddle = complex_exp(t);

	// see if we are dealing with top of bottom of the wing
	// bool is_wingtop = n < N / 2;
	int next_pow = pow(2, stage + 1);
	int cur_pow = pow(2, stage);
	int is_wingtop = 0;
	if (next_pow % n < cur_pow)
	{
		is_wingtop = 1;
	}

	// first stage handled separately
	if (stage == 0)
	{
		// first stage, need to bit reverse indices
		int reversed_n = reversed_bit_indices[n];
		if (is_wingtop == 1)
		{
//			if (n + 1 == N)
//			{
//				m_butterfly_map[stage][n] =
//					{
//					twiddle.real, twiddle.imag, // twiddle
//					reversed_n, reversed_n // save indices
//					};
//				return;
//			}
//			else
			{
				int reverse_n_plus1 = reversed_bit_indices[n + 1];
				m_butterfly_map[stage][n] =
					{
						twiddle.real, twiddle.imag, // twiddle
						reversed_n, reverse_n_plus1 // indices
					};
			}
		}
		else
		{
//			if (n == 0)
//			{
//				m_butterfly_map[stage][n] =
//					{
//					twiddle.real, twiddle.imag, // twiddle
//					reversed_n, reversed_n // save indices
//					};
//			}
//			else
			{
				int reverse_n_minus1 = reversed_bit_indices[n - 1];
				m_butterfly_map[stage][n] =
					{
					twiddle.real, twiddle.imag, // twiddle
					reverse_n_minus1, reversed_n // save indices
					};
			}
		}
	}

	// log2n stages general case
	else
	{
		int butterfly_span = pow(2, stage);

		if (is_wingtop == 1)
		{
			std::cout << "TOP index: " << n + butterfly_span << ", stage: " <<  stage << ", n: " << n << std::endl;
			m_butterfly_map[stage][n] =
				{
				twiddle.real, twiddle.imag, // twiddle
				n, n + butterfly_span // indices
				};
		}
		else
		{
			std::cout << "BOT index: " << n - butterfly_span << ", stage: " <<  stage << ", n: " << n << std::endl;
			m_butterfly_map[stage][n] =
				{
				twiddle.real, twiddle.imag, // twiddle
				n - butterfly_span, n // save indices
				};
		}
	}
}

void fastfft::horizontal_pass
(
	int x,
	int y
)
{
	Eigen::Vector4d butterfly = m_butterfly_map[m_stage][x];
	if (m_pingpong == 0)
	{
		int index1 = butterfly[2];
		int index2 = butterfly[3];
		// TODO, get from pingpong map
		Eigen::Vector2d a = m_pingpong0[index1][y];
		complex p = {a[0], a[1]};
		Eigen::Vector2d b = m_pingpong0[index2][y];
		complex q = {b[0], b[1]};

		complex twiddle = {butterfly[0], butterfly[1]};
		complex m = complex_mul(twiddle, q);
		complex H = complex_add(p, m);

		m_pingpong1[x][y] = {H.real, H.imag};
	}
	else
	{
		int index1 = butterfly[2];
		int index2 = butterfly[3];
		// TODO, get from pingpong map
		Eigen::Vector2d a = m_pingpong1[index1][y];
		complex p = {a[0], a[1]};
		Eigen::Vector2d b = m_pingpong1[index2][y];
		complex q = {b[0], b[1]};

		complex twiddle = {butterfly[0], butterfly[1]};
		complex m = complex_mul(twiddle, q);
		complex H = complex_add(p, m);

		m_pingpong0[x][y] = {H.real, H.imag};
	}
}

void fastfft::vertical_pass
	(
		int x,
		int y
	)
{
	Eigen::Vector4d butterfly = m_butterfly_map[m_stage][y];
	if (m_pingpong == 0)
	{
		int index1 = butterfly[2];
		int index2 = butterfly[3];
		// TODO, get from pingpong map
		Eigen::Vector2d a = m_pingpong0[index1][x];
		complex p = {a[0], a[1]};
		Eigen::Vector2d b = m_pingpong0[index2][x];
		complex q = {b[0], b[1]};

		complex twiddle = {butterfly[0], butterfly[1]};
		complex m = complex_mul(q, twiddle);
		complex H = complex_add(p, m);

		m_pingpong1[x][y] = {H.real, H.imag};
	}
	else
	{
		int index1 = butterfly[2];
		int index2 = butterfly[3];
		// TODO, get from pingpong map
		Eigen::Vector2d a = m_pingpong1[index1][x];
		complex p = {a[0], a[1]};
		Eigen::Vector2d b = m_pingpong1[index2][x];
		complex q = {b[0], b[1]};

		complex twiddle = {butterfly[0], butterfly[1]};
		complex m = complex_mul(q, twiddle);
		complex H = complex_add(p, m);

		m_pingpong0[x][y] = {H.real, H.imag};
	}
}

void fastfft::do_inversion
(
	int x,
	int y
	)
{
	double sign = (x + y) % 2 == 0 ? 1.0 : -1.0;
	// double N_squared = N * N;
	double N_squared = N;

	if (m_pingpong == 0)
	{
		// get the real part of the pingpong
		float real = m_pingpong0[x][y][0];
		float imag = m_pingpong0[x][y][1];
		results[x][y] = {sign * real / N_squared, sign * imag / N_squared};
	}
	else
	{
		// get the real part of the pingpong
		float real = m_pingpong1[x][y][0];
		float imag = m_pingpong1[x][y][1];
		results[x][y] = {sign * real / N_squared, sign * imag / N_squared};
	}
}

std::vector<Eigen::Vector2d> fastfft::do_fft
(
	std::vector<std::vector<Eigen::Vector2d>> input
	)
{
	// set the pingpong texture to the input
	for (int x = 0; x < N; x++)
	{
		for (int y = 0; y < N; y++)
		{
			m_pingpong0[x][y] = input[x][y];
			m_pingpong1[x][y] = input[x][y];
		}
	}
	m_pingpong = 0;

	// do 1D fft for vetical
	for (int stage = 0; stage < log_2_N; stage++)
	{
		m_stage = stage;
		for (int x = 0; x < N; x++)
		{
			for (int y = 0; y < N; y++)
			{
				compute_butterfly(false, x, y);
			}
		}
		m_pingpong = !m_pingpong;
	}

	// do 1D fft for horizontal
	for (int stage = 0; stage < log_2_N; stage++)
	{
		m_stage = stage;
		for (int x = 0; x < N; x++)
		{
			for (int y = 0; y < N; y++)
			{
				compute_butterfly(true, y, x);
			}
		}
		m_pingpong = !m_pingpong;
	}

	// do the inversion for each point
	for (int n = 0; n < N; n++)
	{
		for (int m = 0; m < N; m++)
		{
			do_inversion(n, m);
		}
	}

	// make a 1D vector to return
	std::vector<Eigen::Vector2d> results_1d;
	for (int x = 0; x < N; x++)
	{
		for (int y = 0; y < N; y++)
		{
			results_1d.push_back({results[x][y].real, results[x][y].imag});
		}
	}

	return results_1d;
}

std::vector<Eigen::Vector2d> fastfft::fft2(std::vector<Eigen::Vector2d> input, int stride, int offset)
{
	for (int i = 0; i < N; i++)
	{
		int index = reversed_bit_indices[i] * stride + offset;
		c[which][i] = {input[index][0], input[index][1]};
	}

	int loops       = N>>1;
	int size        = 1<<1;
	int size_over_2 = 1;
	int w_          = 0;
	for (int i = 1; i <= log_2_N; i++) {
		which ^= 1;
		for (int j = 0; j < loops; j++) {
			for (int k = 0; k < size_over_2; k++) {
				c[which][size * j + k] =
					complex_add(c[which][size * j + k],
					complex_mul(c[which^1][size * j + size_over_2 + k], m_butterfly[w_][k]));
			}

			for (int k = size_over_2; k < size; k++) {
				c[which][size * j + k] =
					complex_sub(c[which^1][size * j - size_over_2 + k],
					complex_mul(c[which^1][size * j + k], m_butterfly[w_][k - size_over_2]));
			}
		}
		loops       >>= 1;
		size        <<= 1;
		size_over_2 <<= 1;
		w_++;
	}

	std::vector<Eigen::Vector2d> output = input;
	for (int i = 0; i < N; i++)
	{
		output[i * stride + offset] = { c[which][i].real, c[which][i].imag };
	}
	return output;
}

