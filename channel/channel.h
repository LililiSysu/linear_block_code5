/*****************************************************************//**
 * \file   channel.h
 * \brief  including class of AWGN, BPSK, GF_trans, Measure
 * 
 * \author 26259
 * \date   October 2023
 *********************************************************************/

#pragma once

#include"../my_lib/Matrix.h"
#include"../GF/GF2.h"
#include<random>

// we use my_double only in channel modeling, be careful
#ifdef use_my_double
#include"../my_lib/my_double.h"
#else
#define my_double double
#endif

// do not know why we cannot seperate .h file and .cpp file

class AWGN {

public:
	static my_double sigma;

	/**
	 * .pass this AWGN channel
	 *
	 * \param c: vector to pass this channel
	 * \param _sigma: if given, reset sigma to generate noise, else use old sigma
	 * \return the vector after this channel
	 */
	static Matrix<my_double> pass(const Matrix<my_double>& c, my_double _sigma = -1) {
		sigma = _sigma == -1 ? sigma : _sigma;			// reseting sigma if given _sigma

		Matrix<my_double> result(c.row(), c.col());
		int result_size = result.size();
		for (int i = 0; i < result_size; ++i) {
			result(i) = (double)c(i) + my::rand_ga() * (double)sigma;		// prevent counting double operation
		}
		return result;
	}
	
	/**
	 * .pass this AWGN channel, generate the Gaussian noise by standard random function
	 * it seems no difference, over distribution, with our own written gaussian random varible
	 * when number of random variable generated large enough, but this funciton is more precise.
	 * please use this during simulation of AWGN channel. Using this standard gaussian generator
	 * will include higher error rate than our own written gaussian generator. be careful.
	 * 
	 * it will be normal that the error rate induce by this function is higher than that 
	 * published over the internet, please trust this function anyway, and don't be surprise.
	 * 
	 * the time consume of this funciton is alomst the same as 'pass'
	 * 
	 * \param c: vector to pass this channel
	 * \param _sigma: if given, reset sigma to generate noise, else use old sigma
	 * \return the vector after this channel
	 */
	static Matrix<my_double> pass_standard(const Matrix<my_double>& c, my_double _sigma = -1) {
		sigma = _sigma == -1 ? sigma : _sigma;			// reseting sigma if given _sigma

		random_device rd;
		mt19937 gen(rd());
		normal_distribution<double> normal(0.0, (double)sigma);

		Matrix<my_double> result(c.row(), c.col());
		int result_size = result.size();
		for (int i = 0; i < result_size; ++i) {
			result(i) = (double)c(i) + normal(gen);		// prevent counting double opeartion
		}
		return result;
	}
};
my_double AWGN::sigma = 1;// set sigma be 1 as default, statandard Gaussian distribution

class BPSK {
public:
	/**
	 * .perform BPSK modulation, 0 to 1 and 1 to -1
	 *
	 * \param v: vector to modulate
	 * \return BPSK modulated vector
	 */
	 static Matrix<my_double> modulation(const Matrix<GF2>& v) {
		Matrix<my_double> result(v.row(), v.col());
		for (int i = 0; i < result.size(); ++i) {
			result(i) = v(i) == 0 ? 1 : -1;
		}
		return result;
	}

	/**
	 * .perform BPSK demodulation, to 0 if symbol > 0, else to 1
	 *
	 * \param v: vector to demodulate
	 * \return BPSK demodulated vector
	 */
	static Matrix<GF2> demodulation(const Matrix<my_double>& r) {
		Matrix<GF2> result(r.row(), r.col());
		for (int i = 0; i < result.size(); ++i) {
			result(i) = r(i) > 0 ? 0 : 1;
		}
		return result;
	}

	static my_double modulation(const GF2& v) {
		return v == 0 ? 1 : -1;
	}
};

template<class T, int bits_per_symbol>
class GF_trans {
public:
	/**
	 * .turn a vector on field T to GF2
	 */
	Matrix<GF2> to_bits(const Matrix<T>& v) {
		int sss = v.size();
		Matrix<GF2> result(1, sss * bits_per_symbol);
		int num_symbols = 0;
		for (int i = 0; i < sss; ++i) {
			int sym = (int)v(i);
			for (int j = 0; j < bits_per_symbol; ++j) {
				result(num_symbols + j) = sym & 1;
				sym >>= 1;
			}
			num_symbols += bits_per_symbol;
		}
		return result;
	}

	Matrix<T> to_symbol(const Matrix<GF2>& x) {
		int sss = x.size();
		int symbol_sss = sss / bits_per_symbol;
		Matrix<T> result(1, symbol_sss);
		int symbol_pos = 0;
		for (int i = 0; i < symbol_sss; ++i) {
			int sym = 0;
			for (int j = bits_per_symbol -1; j >=0; --j) {
				sym <<= 1;
				sym = sym | (int)x(symbol_pos + j);
			}
			result(i) = T(sym);
			symbol_pos += bits_per_symbol;
		}
		return result;
	}
};

class Measure {
public:
	static int error_num;

	/**
	 * .
	 * 
	 * \param r: soft sequence
	 * \param v: hard sequence, soft to hard mapping relation is, if soft_val >= 0, hard_val = 0  else hard_val = 1
	 * \param error_num: the number of different symbols between hard decoded r and v
	 * \return correlation disvrepancy
	 */
	static my_double correlation_discrepancy_v(const Matrix<my_double>& r, const Matrix<GF2>& v) {
		my_double result = 0;
		error_num = 0;
		int len = r.size();

		for (int i = 0; i < len; ++i) {
			if ((r(i) >= 0 && v(i) == 0) || (r(i) < 0 && v(i) == 1));
			else {
				result += my::abs(r(i));
				error_num++;
			}
		}
		return result;
	}

	/**
	 * .this for faster computation
	 *
	 * \param r_abs: abs value soft sequence
	 * \param hdr: the hard decoded soft sequence, under BPSK demodulation
	 * \param v: hard sequence, soft to hard mapping relation is, if soft_val >= 0, hard_val = 0  else hard_val = 1
	 * \param error_num: the number of different symbols between hdr and v
	 * \return correlation disvrepancy
	 */
	static my_double correlation_discrepancy_v(const Matrix<my_double>& r_abs, const Matrix<GF2>& hdr, const Matrix<GF2>& v) {

		// not important for speed
		my_double result = 0;
		error_num = 0;
		int len = r_abs.size();

		for (int i = 0; i < len; ++i) {
			result += hdr(i) != v(i) ? r_abs(i) : 0;
			//error_num += hdr(i) != v(i) ? 1 : 0;
			error_num += hdr(i) != v(i);
		}
		return result;
	}

	/**
	 * .Euclidean distance between 2 real vector 'r' and 'v', also accumulating error_num as |r(i)-v(i)| > 1
	 */
	static my_double Euclidean_distance(const Matrix<my_double>& r, const Matrix<my_double>& v) {
		my_double result = 0;
		int len = r.size();
		my_double tmp; 
		error_num = 0;

		for (int i = 0; i < len; ++i) {
			tmp = r(i) - v(i);
			result += tmp * tmp;
			error_num += my::abs(tmp) > 1;
		}
		return result;
	}

	/**
	 * .Euclidean distance between 'r' and BPSK::modulation('v')
	 * 
	 * \param r: real vector
	 * \param v: GF2 vector
	 * \param n_plus_sum_r_squared_minus_2_r_abs: n + sum of (r(i)^2-2|r(i)|)
	 * \return 
	 */
	static my_double Euclidean_distance(const Matrix<my_double>& r, const Matrix<GF2>& v, my_double n_plus_sum_r_squared_minus_2_r_abs) {
		return n_plus_sum_r_squared_minus_2_r_abs + 4 * correlation_discrepancy_v(r, v);
	}

	/**
	 * .
	 *
	 * \param r: soft sequence
	 * \param v: hard sequence, soft to hard mapping relation is, if soft_val >= 0, hard_val = 0  else hard_val = 1
	 * \param error_num: the number of different symbols between hard decoded r and v
	 * \return correlation disvrepancy
	 */
	static double correlation_discrepancy_v_no_count_ope(const Matrix<my_double>& r, const Matrix<GF2>& v) {
		double result = 0;
		error_num = 0;
		int len = r.size();

		for (int i = 0; i < len; ++i) {
			if (((double)r(i) >= 0 && (int)v(i) == 0) || ((double)r(i) < 0 && (int)v(i) == 1));
			else {
				result += my::abs((double) r(i));
				error_num++;
			}
		}
		return result;
	}
};

int Measure::error_num = 0;
