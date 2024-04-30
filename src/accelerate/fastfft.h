//
// Created by Michael Foiani on 4/29/24.
//

#ifndef FASTFFT_H
#define FASTFFT_H

#include <vector>
#include "Eigen/Dense"

struct complex
{
	double real;
	double imag;
};

class fastfft
{
public:
	fastfft(int N);

	// INSTANCE VARIABLES
	// Ocean grid params
	// unsigned int N; // total number of grid points in one D! must be power of 2 and square root must be an integer
	// int log2_N; // log base 2 of N

	unsigned int N, which, log_2_N;
	float pi2;
	unsigned int *reversed;
	std::vector<std::vector<complex>> c;


	// Final results from the computation
	std::vector<std::vector<complex>> results;
	std::vector<Eigen::Vector2d> do_fft(
		std::vector<std::vector<Eigen::Vector2d>> input
	);

	// underlying wavevectors
	std::vector<Eigen::Vector2d> h0_k;
	std::vector<Eigen::Vector2d> h0_minsk;

	void do_inversion
		(
			int x,
			int y
		);

	// fft fast butterfly computations
	void vertical_pass
		(
			int x,
			int y
		);
	void horizontal_pass
		(
			int x,
			int y
		);
	inline void compute_butterfly(bool is_vertical, int x, int y)
	{
		if (is_vertical)
		{
			vertical_pass(x, y);
		}
		else
		{
			horizontal_pass(x, y);
		}
	}

	std::vector<Eigen::Vector2d> fft2(std::vector<Eigen::Vector2d> input, int stride, int offset);

	// COMPLEX OPERATIONS
	static inline complex complex_add(complex a, complex b)
	{
		return {a.real + b.real, a.imag + b.imag};
	}
	static inline complex complex_sub(complex a, complex b)
	{
		return {a.real - b.real, a.imag - b.imag};
	}
	static inline complex complex_mul(complex a, complex b)
	{
		return {a.real * b.real - a.imag * b.imag, a.real * b.imag + a.imag * b.real};
	}
	static inline complex complex_exp(double x)
	{
		return {cos(x), sin(x)};
	}
	static inline complex complex_conj(complex x)
	{
		return {x.real, -x.imag};
	}

private:
	// butterfly map for fast fft
	std::vector<std::vector<Eigen::Vector4d>> m_butterfly_map;
	std::vector<std::vector<complex>> m_butterfly;
	std::vector<int> reversed_bit_indices;
	void precompute_butterfly_map(
		int log2_N,
		int sqrt_N
		);

	// pingpong texture for the 1D decompoisitons
	int m_pingpong; // current resulting texture in use
	std::vector<std::vector<Eigen::Vector2d>> m_pingpong0;
	std::vector<std::vector<Eigen::Vector2d>> m_pingpong1;
	int m_stage = 0; // current stage of the fft compositions

};


#endif //FASTFFT_H
