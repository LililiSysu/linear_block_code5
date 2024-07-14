/*****************************************************************//**
 * \file   simulation.h
 * \brief  simulation structure for OSD_Chase_orig.h in 'test' class
 * 
 * \author 26259
 * \date   October 2023
 *********************************************************************/

#pragma once

#include"../my_lib/Matrix.h"
#include"../GF/GF2e.h"
#include"../code/BCH.h"
#include"../code/eBCH.h"
#include"../code/RS.h"
#include"../code/Weight_Spectrum.h"
#include"channel.h"
#include"../reprocess/Hybrid_Chase2_OSD.h"
#include<ctime>

// we use my_double only in simulation, be careful
#ifdef use_my_double
#include"../my_lib/my_double.h"
#else
#define my_double double
#endif  // use_my_double

using namespace std;

enum class decode_type{
	Algebra,
	GMD,
	Chase2,
	Chase3,
	OSD,
	Hybrid_Chase2_OSD,
};

enum class RS_decode_type {
	BM,
	GS
};

struct simulation_param {
	decode_type dt;
	int param_1;
	int param_2;

	// for Hybrid_Chase2_OSD decoding, param_1 is OSD order, param_2 is Chase2 flip bit

	simulation_param(decode_type _dt = decode_type::Algebra, int _param_1 = 0, int _param_2 = 0) {
		dt = _dt;
		param_1 = _param_1;
		param_2 = _param_2;
	}
};

class simulation
{
public:
	static int _iteration;
	static int _error_frame_upper_bound;

	static my_double BER;
	static my_double FER;
	static int frame_proceed;

	static my_double GF2e_operation_cnt_ave;		// average GF2e operation during decoding
	static my_double GF2_operation_cnt_ave;			// average GF2 operation during decoding
	static my_double my_double_operation_cnt_ave;	// average my_double operation during decoding
	static double time_consume_ave;					// average time consume

	/**
	 * .simulate decoding to get BER-SNR curve line
	 * 
	 * \param code: type of code, such as BCH<4,2>
	 * \param SNR_dB: SNR=Eb/N0 in dB, energy per information bit over noise power specturm
	 * \param iteration: number of encode and decode for simulation
	 * \param dt: decode type, can be one of 
	 *		decode_type::Algebra, decode_type::GMD, decode_type::Chase2, decode_type::Chase3, decode_type::OSD
	 * \param dt_para: decode type parameter, if choosen decode_type::Chase2, 'dt_para' will be number of bit to flip, 
	 * if choosen decode_type::OSD, 'dt_para' will be order of OSD
	 * 
	 * \return  please fetch number from static varible of this class, including FER(Frame Error Rate), BER(Bit Error,Rate), 
	 * frame_proceed, GF2e_operation_cnt_ave and GF2_operation_cnt_ave
	 */	
	template<class code_type> 
	static void	linear_block_code(code_type& code, my_double SNR_dB, int iteration = 1, int error_frame_upper_bound = 100, \
			decode_type dt = decode_type::Algebra, int dt_para_eta = -1, int dt_para_order = -1) {

		//const int error_frame_upper_bound = _error_frame_upper_bound;
		int k = code.get_k();
		int n = code.get_n();
		my_double rate = (my_double)k / n;

		Matrix<GF2> u(1, k, 'b');
		Matrix<GF2> v = code.encode(u);
		Matrix<my_double> c = BPSK::modulation(v);
		Matrix<my_double> e(1, n);

		//cout << "SNR_dB=" << SNR_dB << endl;
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * rate * SNR));	// BPSK modulation {-1,1} has energy 1

		//cout << "sigma=" << sigma << endl;
		int error_symbols = 0;
		int error_Frame = 0; 

		// reset varible
		GF2e_operation_cnt_ave = 0;
		GF2_operation_cnt_ave = 0;

		for (int iter = 0; iter < iteration; ++iter) {
			//cout << "iteration=" << iter << endl;
			for (int i = 0; i < n; ++i) {
				e(i) = my::rand_ga() * sigma;
			}
			//cout << "e" << e;
			// to do: reduce space declearation and deletion, avoid matrix created inside for loop, speed up

			//Matrix<my_double> r = AWGN::pass(c);// with standard Gaussian distribution as default, N(0,1)
			Matrix<my_double> r = e + c;
			//Matrix<my_double> r =c;

			Matrix<GF2> hdr = BPSK::demodulation(r);		// consider use of reference induction
			//cout << "hdr=" << hdr;
			Matrix<GF2> du;

#ifdef count_operation_number
			// reset to zero to prevent overflow
			GF2e_auxiliary_storage::operation_number = 0;
			GF2_auxiliary_storage::operation_number = 0;
			my_double_auxiliary_storage::operation_number = 0;

			unsigned long long GF2e_operation_cnt_before_decode = GF2e_auxiliary_storage::operation_number;
			unsigned long long GF2_operation_cnt_before_decode = GF2_auxiliary_storage::operation_number;
			unsigned long long my_double_operation_cnt_before_decode = my_double_auxiliary_storage::operation_number;

			clock_t start_sub, end_sub;
			start_sub = clock();
#endif // count_operation_number

			switch (dt) {
				case(decode_type::Algebra):du = code.decode(hdr); break;
				case(decode_type::GMD):du = GMD::decode(code, r); break;
				case(decode_type::Chase2):du = Chase::decode2(code, r, dt_para_eta); break;	// flip floor of dmin/2 LRP as default
				case(decode_type::Chase3): du = Chase::decode3(code, r); break;					// worse than GMD::decode
				case(decode_type::OSD):du = OSD::decode(code, r, dt_para_order); break;			// order = dmin/4 as the default
				case(decode_type::Hybrid_Chase2_OSD):du = Hybrid_Chase2_OSD::decode_v(code, r, dt_para_eta, dt_para_order); break;
					// flip dt_para and order dt_para
				default: cout << "decode type input error" << endl;
			}


#ifdef count_operation_number
			end_sub = clock();
			time_consume_ave += ((double)end_sub - start_sub) / CLOCKS_PER_SEC;

			GF2e_operation_cnt_ave += (double) GF2e_auxiliary_storage::operation_number - GF2e_operation_cnt_before_decode;
			GF2_operation_cnt_ave += (double) GF2_auxiliary_storage::operation_number - GF2_operation_cnt_before_decode;
			my_double_operation_cnt_ave += (double) my_double_auxiliary_storage::operation_number - my_double_operation_cnt_before_decode;
#endif // count_operation_number

			/*if (iter % 2000 != 0);
			else{
				if (ret_BER)
					cout << "iter=" << iter << ", bit error rate=" << (my_double)error_symbols / ((my_double)(iter + 1.0) * k) << endl;
				else
					cout << "iter=" << iter << ", frame error rate=" << (my_double)error_Frame / (my_double)(iter + 1.0) << endl;
			}*/

			/* compute error element */
			int err_tmp = u.Hamming_distance(du);
			error_symbols += err_tmp;
			error_Frame += err_tmp != 0;
			if (error_Frame <= error_frame_upper_bound);		// if error_frame_upper_bound frames occur error, return
			else{
				//cout << "(" << SNR_dB << " dB) real iteration=" << iter + 1 << endl;
				BER=(my_double)error_symbols / ((my_double)(iter + 1.0) * k);		// return bit error rate
				FER=(my_double)error_Frame / (my_double)(iter + 1.0);
				frame_proceed = iter + 1;
				GF2e_operation_cnt_ave /= frame_proceed;
				GF2_operation_cnt_ave /= frame_proceed;
				my_double_operation_cnt_ave /= frame_proceed;
				time_consume_ave /= frame_proceed;
				return;
			}
		}

		//cout << "(" << SNR_dB << " dB) real iteration=" << iteration << endl;
		BER = (my_double)error_symbols / ((my_double)(iteration) * k);		// return bit error rate
		FER = (my_double)error_Frame / (my_double)(iteration);
		frame_proceed = iteration;

		GF2e_operation_cnt_ave /= frame_proceed;
		GF2_operation_cnt_ave /= frame_proceed;
		my_double_operation_cnt_ave /= frame_proceed;
		time_consume_ave /= frame_proceed;

		//if (ret_BER)		// throw away
		//	return make_pair((my_double)error_symbols / ((my_double)iteration * k), iteration);		// return bit error rate
		//else
		//	return make_pair((my_double)error_Frame / (my_double)iteration, iteration);
	}

	// warping of linear_block_code
	template<class code_type> 
	static void	linear_block_code_param_struct(code_type& code, simulation_param sm, my_double SNR_dB) {

		switch (sm.dt) {
			// Notice that in 'linear_block_code' last parameter is OSD order and the second last is Chase2 flip bit
		case(decode_type::Algebra):linear_block_code(code, SNR_dB, _iteration, _error_frame_upper_bound, sm.dt, 0, 0); break;
		case(decode_type::GMD):linear_block_code(code, SNR_dB, _iteration, _error_frame_upper_bound, sm.dt, 0, 0); break;
		case(decode_type::Chase3):linear_block_code(code, SNR_dB, _iteration, _error_frame_upper_bound, sm.dt, 0, 0); break;
		case(decode_type::Chase2):linear_block_code(code, SNR_dB, _iteration, _error_frame_upper_bound, sm.dt, sm.param_1, 0); break;
		case(decode_type::OSD):linear_block_code(code, SNR_dB, _iteration, _error_frame_upper_bound, sm.dt, 0, sm.param_1); break;
				// for case of OSD, be careful
		case(decode_type::Hybrid_Chase2_OSD):
			linear_block_code(code, SNR_dB, _iteration, _error_frame_upper_bound, sm.dt, sm.param_2, sm.param_1);
			break;

		default: cout << "decode type input error" << endl;
		}
	}


	/**
	 * .simulate decoding to get BER-SNR curve line, special for RS code
	 *
	 * \param SNR_dB: SNR=Eb/N0 in dB, energy per information bit over noise power specturm
	 * \param dt: decode type, can be one of RS_decode_type::BM, RS_decode_type::GS
	 * \param multiplicity: decode type parameter, if choosen RS_decode_type::GS, 'dt_para' will be multiplicity 
	 *
	 * \return  please fetch number from static varible of this class, including FER(Frame Error Rate), BER(Bit Error,Rate),
	 * frame_proceed, GF2e_operation_cnt_ave and GF2_operation_cnt_ave, my_double_operation_cnt_ave, time_consume_ave
	 */
	template<int m, int k>
	static void	RS_code(my_double SNR_dB, RS_decode_type dt = RS_decode_type::BM, int multiplicity = 1) {

		RS<m, k> code;
		GF_trans<GF2e<m>, m> gt;

		int n = code.get_n();
		my_double rate = (my_double)k / n;

		Matrix<GF2e<m>> u(1, k, '1');		// all one symbol transmitted
		Matrix<GF2e<m>> v_sym = code.encode_by_code_locators(u);
		Matrix<GF2> v = gt.to_bits(v_sym);
		Matrix<my_double> c = BPSK::modulation(v);
		int bit_num = c.size();
		//cout << "bit_num = " << bit_num << endl;
		Matrix<my_double> e(1, bit_num);

		//cout << "SNR_dB=" << SNR_dB << endl;
		my_double SNR = pow(10, SNR_dB / 10);
		//cout << "SNR=" << SNR << endl;
		my_double sigma = sqrt(1.0 / (2 * rate * SNR));	// BPSK modulation {-1,1} has energy 1

		//cout << "sigma=" << sigma << endl;
		int error_symbols = 0;
		int error_Frame = 0;

		// reset varible
		GF2e_operation_cnt_ave = 0;
		GF2_operation_cnt_ave = 0;

		for (int iter = 0; iter < _iteration; ++iter) {
			//cout << "_iteration=" << iter << endl;
			for (int i = 0; i < bit_num; ++i) {
				e(i) = my::rand_ga() * sigma;
			}
			//cout << "e" << e;
			// to do: reduce space declearation and deletion, avoid matrix created inside for loop, speed up

			//Matrix<my_double> r = AWGN::pass(c);// with standard Gaussian distribution as default, N(0,1)
			Matrix<my_double> r = e + c;
			//Matrix<my_double> r =c;

			Matrix<GF2> hdr = BPSK::demodulation(r);		// consider use of reference induction

			//cout << "hdr=" << hdr;
			Matrix<GF2e<m>> hdr_sym = gt.to_symbol(hdr);
			Matrix<GF2e<m>> du;

#ifdef count_operation_number
			// reset to zero to prevent overflow
			GF2e_auxiliary_storage::operation_number = 0;
			GF2_auxiliary_storage::operation_number = 0;
			my_double_auxiliary_storage::operation_number = 0;

			unsigned long long GF2e_operation_cnt_before_decode = GF2e_auxiliary_storage::operation_number;
			unsigned long long GF2_operation_cnt_before_decode = GF2_auxiliary_storage::operation_number;
			unsigned long long my_double_operation_cnt_before_decode = my_double_auxiliary_storage::operation_number;

			clock_t start_sub, end_sub;
			start_sub = clock();
#endif // count_operation_number

			switch (dt) {
			case(RS_decode_type::BM):du = code.decode_BM_by_code_locators(hdr_sym); break;
			case(RS_decode_type::GS):du = code.decode_GS_by_code_locators(hdr_sym, multiplicity); break;
			default: cout << "decode type input error" << endl;
			}


#ifdef count_operation_number
			end_sub = clock();
			time_consume_ave += ((double)end_sub - start_sub) / CLOCKS_PER_SEC;

			GF2e_operation_cnt_ave += (double) GF2e_auxiliary_storage::operation_number - GF2e_operation_cnt_before_decode;
			GF2_operation_cnt_ave += (double) GF2_auxiliary_storage::operation_number - GF2_operation_cnt_before_decode;
			my_double_operation_cnt_ave += (double) my_double_auxiliary_storage::operation_number - my_double_operation_cnt_before_decode;
#endif // count_operation_number

			/*if (iter % 2000 != 0);
			else{
				if (ret_BER)
					cout << "iter=" << iter << ", bit error rate=" << (my_double)error_symbols / ((my_double)(iter + 1.0) * k) << endl;
				else
					cout << "iter=" << iter << ", frame error rate=" << (my_double)error_Frame / (my_double)(iter + 1.0) << endl;
			}*/

			/* compute error element */
			
			//cout << "du" << du;
			int err_tmp = u.Hamming_distance(du);
			error_symbols += err_tmp;
			error_Frame += err_tmp != 0;
			if (error_Frame <= _error_frame_upper_bound);		// if _error_frame_upper_bound frames occur error, return
			else {
				//cout << "(" << SNR_dB << " dB) real _iteration=" << iter + 1 << endl;
				BER = (my_double)error_symbols / ((my_double)(iter + 1.0) * k);		// return bit error rate
				FER = (my_double)error_Frame / (my_double)(iter + 1.0);
				frame_proceed = iter + 1;
				GF2e_operation_cnt_ave /= frame_proceed;
				GF2_operation_cnt_ave /= frame_proceed;
				my_double_operation_cnt_ave /= frame_proceed;
				time_consume_ave /= frame_proceed;
				return;
			}
		}

		//cout << "(" << SNR_dB << " dB) real _iteration=" << _iteration << endl;
		BER = (my_double)error_symbols / ((my_double)(_iteration)*k);		// return bit error rate
		FER = (my_double)error_Frame / (my_double)(_iteration);
		frame_proceed = _iteration;

		GF2e_operation_cnt_ave /= frame_proceed;
		GF2_operation_cnt_ave /= frame_proceed;
		my_double_operation_cnt_ave /= frame_proceed;
		time_consume_ave /= frame_proceed;

		//if (ret_BER)		// throw away
		//	return make_pair((my_double)error_symbols / ((my_double)_iteration * k), _iteration);		// return bit error rate
		//else
		//	return make_pair((my_double)error_Frame / (my_double)_iteration, _iteration);
	}
};

int simulation::_iteration = 5000000;		// change to 600 if test operation number
int simulation::_error_frame_upper_bound = 100;

my_double simulation::FER = 0;
my_double simulation::BER = 0;
int simulation::frame_proceed = 0;

my_double simulation::GF2_operation_cnt_ave = 0;
my_double simulation::GF2e_operation_cnt_ave = 0;
my_double simulation::my_double_operation_cnt_ave = 0;
double simulation::time_consume_ave = 0;
