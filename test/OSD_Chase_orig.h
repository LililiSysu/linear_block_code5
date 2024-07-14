/*****************************************************************//**
 * \file   OSD_Chase_orig.h
 * \brief  simulation for OSD and Chase algorithm with an original structure. No use.now
 * 
 * \author 26259
 * \date   October 2023
 *********************************************************************/

#pragma once

#include"test_common.h"

class test_OSD_Chase_orig {
private:

	template<class code_type>
	static void simulate_param(code_type& code, simulation_param sm, const Matrix<my_double>& SNR_db) {

		// print simulation info
		{
			code.print_info();
			switch (sm.dt) {
			case(decode_type::Algebra):cout << "Algebra simulation" << endl; break;
			case(decode_type::GMD):cout << "GMD simulation" << endl; break;
			case(decode_type::Chase2):cout << "Chase2 simulation, flip bit = " << sm.param_1 << endl; break;
			case(decode_type::Chase3): cout << "Chase3 simulation" << endl; break;
			case(decode_type::OSD):cout << "OSD simulation, order = " << sm.param_1 << endl; break;
			case(decode_type::Hybrid_Chase2_OSD):
				cout << "Hybrid_Chase2_OSD simulation, order = " << sm.param_1 << ", eta = " << sm.param_2 << endl; break;
			default: cout << "decode type input error" << endl;
			}

			cout << "iteration frames = " << simulation::_iteration << endl;
			cout << "error frames = " << simulation::_error_frame_upper_bound << endl;
		}

		// simulation result
		const int SNR_num = SNR_db.size();

		Matrix<my_double> BER(1, SNR_num);
		Matrix<my_double> FER(1, SNR_num);
		Matrix<int> Frame(1, SNR_num);
		Matrix<my_double> GF2_operation(1, SNR_num);
		Matrix<my_double> GF2e_operation(1, SNR_num);
		Matrix<my_double> my_double_operation(1, SNR_num);
		Matrix<my_double> time_consume_ave(1, SNR_num);

		clock_t start, end;
		start = clock();		// starting time

		for (int i = 0; i < SNR_num; i++) {
			simulation::linear_block_code_param_struct(code, sm, SNR_db(i));
			BER(i) = simulation::BER;
			FER(i) = simulation::FER;
			Frame(i) = simulation::frame_proceed;
			GF2_operation(i) = simulation::GF2_operation_cnt_ave;
			GF2e_operation(i) = simulation::GF2e_operation_cnt_ave;
			my_double_operation(i) = simulation::my_double_operation_cnt_ave;
			time_consume_ave(i) = simulation::time_consume_ave;
		}

		end = clock();

		// print simulation result
		{
			cout << "\nSNR_db" << SNR_db;

#ifdef count_operation_number
			cout << "GF2_operation" << GF2_operation;
			cout << "GF2e_operation" << GF2e_operation;
			cout << "my_double_operation" << my_double_operation;

			cout.setf(ios::scientific);
			cout << "time_consume_ave" << time_consume_ave;
			cout.unsetf(ios::floatfield);
#else
			cout.setf(ios::scientific);
			cout << "BER" << BER;
			cout << "FER" << FER;
			cout.unsetf(ios::floatfield);

			cout << "Frame" << Frame;
#endif

			cout << "\ntotal time consume: " << ((double)end - start) / CLOCKS_PER_SEC << "s\n" \
				<< "\n------------------------------\n------------------------------\n\n";
		}
	}

	template<class code_type>
	static void simulate_param_auto(code_type& code, simulation_param sm, my_double SNR_start = 0.0, my_double SNR_step = 0.5) {
		// print code info
		{
			code.print_info();
			switch (sm.dt) {
			case(decode_type::Algebra):cout << "Algebra simulation" << endl; break;
			case(decode_type::GMD):cout << "GMD simulation" << endl; break;
			case(decode_type::Chase2):cout << "Chase2 simulation, flip bit = " << sm.param_1 << endl; break;
			case(decode_type::Chase3): cout << "Chase3 simulation" << endl; break;
			case(decode_type::OSD):cout << "OSD simulation, order = " << sm.param_1 << endl; break;
			case(decode_type::Hybrid_Chase2_OSD):
				cout << "Hybrid_Chase2_OSD simulation, order = " << sm.param_1 << ", eta = " << sm.param_2 << endl; break;
			default: cout << "decode type input error" << endl;
			}

			cout << "iteration frames = " << simulation::_iteration << endl;
			cout << "error frames = " << simulation::_error_frame_upper_bound << endl;
		}

		const int SNR_num = 30;			// large enougth SNR number
		Matrix<my_double> SNR_db(1, SNR_num, 'v');

		Matrix<my_double> BER(1, SNR_num, 'v');
		Matrix<my_double> FER(1, SNR_num, 'v');
		Matrix<int> Frame(1, SNR_num, 'v');
		Matrix<my_double> GF2_operation(1, SNR_num, 'v');
		Matrix<my_double> GF2e_operation(1, SNR_num, 'v');
		Matrix<my_double> my_double_operation(1, SNR_num, 'v');
		Matrix<my_double> time_consume_ave(1, SNR_num, 'v');		// simulation result

		clock_t start, end;
		start = clock();		// starting time

		SNR_db.push_back(SNR_start);
		do {
			simulation::linear_block_code_param_struct(code, sm, SNR_db.back());
			BER.push_back(simulation::BER);
			FER.push_back(simulation::FER);
			Frame.push_back(simulation::frame_proceed);
			GF2_operation.push_back(simulation::GF2_operation_cnt_ave);
			GF2e_operation.push_back(simulation::GF2e_operation_cnt_ave);
			my_double_operation.push_back(simulation::my_double_operation_cnt_ave);
			time_consume_ave.push_back(simulation::time_consume_ave);

			SNR_db.push_back(SNR_db.back() + SNR_step);
		} while (Frame.back() != simulation::_iteration);
		SNR_db.pop_back();

		end = clock();

		// print simulation result
		{
			cout << "\nSNR_db" << SNR_db;

#ifdef count_operation_number
			cout << "GF2_operation" << GF2_operation;
			cout << "GF2e_operation" << GF2e_operation;
			cout << "my_double_operation" << my_double_operation;

			cout.setf(ios::scientific);
			cout << "time_consume_ave" << time_consume_ave;
			cout.unsetf(ios::floatfield);
#else
			cout.setf(ios::scientific);
			cout << "BER" << BER;
			cout << "FER" << FER;
			cout.unsetf(ios::floatfield);

			cout << "Frame" << Frame;
#endif

			cout << "\ntotal time consume: " << ((double)end - start) / CLOCKS_PER_SEC << "s\n" \
				<< "\n------------------------------\n------------------------------\n\n";
		}
	}

	template<class code_type>
	static void simulate_operation(code_type& code) {

		// set cout into file

		{
			// print all to file
			char file_name[55] = { 0 };
#ifdef count_operation_number
			sprintf_s(file_name, 55, "%d-%d-operation.txt", code.get_n(), code.get_k());
#else
			sprintf_s(file_name, 55, "%d-%d-FER.txt", code.get_n(), code.get_k());
#endif
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}/**/


		Matrix<my_double> SNR_db(0, 0.5, 8, 'd');
		simulation::_iteration = 1000;
		simulation::_error_frame_upper_bound = 1000;

		for (int i = 0; i < (int)simulation_param_set_general.size(); ++i) {
			simulate_param(code, simulation_param_set_general[i], SNR_db);
		}
		fclose(stdout);
	}

	template<class code_type>
	static void simulate_FER_auto(code_type& code, my_double SNR_start = 0.0, my_double SNR_step = 0.5) {
		// set cout into file
		{
			// print all to file
			char file_name[55] = { 0 };
#ifdef count_operation_number
			sprintf_s(file_name, 55, "%d-%d-operation.txt", code.get_n(), code.get_k());
#else
			sprintf_s(file_name, 55, "%d-%d-FER.txt", code.get_n(), code.get_k());
#endif
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		simulation::_iteration = 5000000;
		simulation::_error_frame_upper_bound = 100;

		for (int i = 0; i < (int)simulation_param_set_general.size(); ++i) {
			simulate_param_auto(code, simulation_param_set_general[i], SNR_start, SNR_step);
		}
		fclose(stdout);
	}

	template<class code_type>
	static void simulate_FER_specified_SNR(code_type& code, const Matrix<my_double>& SNR_db) {

		// set cout into file
		{
			// print all to file
			char file_name[55] = { 0 };
#ifdef count_operation_number
			sprintf_s(file_name, 55, "%d-%d-operation.txt", code.get_n(), code.get_k());
#else
			sprintf_s(file_name, 55, "%d-%d-FER.txt", code.get_n(), code.get_k());
#endif
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}

		simulation::_iteration = 5000000;
		simulation::_error_frame_upper_bound = 100;

		for (int i = 0; i < (int)simulation_param_set_general.size(); ++i) {
			simulate_param(code, simulation_param_set_general[i], SNR_db);
		}
		fclose(stdout);
	}

public:

	static const int _m_ = 7;
	static const int _t_ = 10;

	static const int _m1_ = 6;
	static const int _t11_ = 5;
	static const int _t12_ = 4;
	static const int _t13_ = 3;

	static const int _m2_ = 7;
	static const int _t21_ = 6;

	static const int _m3_ = 8;
	static const int _t31_ = 11;
	static const int _t32_ = 5;
	static const int _t33_ = 4;

	/* matrix class */
	static void matrix_identity_test() {
		Matrix<int> A(4, 4, 'i');
		cout << "A" << A;
	}
	static void matrix_mul() {

		typedef my_double Ty;
		//Ty::init();
		cout << is_pod<int>::value << endl;
		cout << is_pod<Ty>::value << endl;	// return 1, its declearation needn't construction fucntioin and is faster
		clock_t start, end;
		const int n = 1000;
		int n2 = n * n;
		Matrix<Ty> m1(n, n);
		for (int i = 0; i < n2; ++i) {
			m1(i) = my::rand_u();
		}
		Matrix<Ty> m2(n, n);
		for (int i = 0; i < n2; ++i) {
			m2(i) = my::rand_u();
		}
		Matrix<Ty> m3(n, n);				// for my_double 
		for (int i = 0; i < n2; ++i) {
			m3(i) = my::rand_u();
		}

		/*Matrix<Ty> m1(dim, dim, 'b');
		Matrix<Ty> m2(dim, dim, 'b');
		Matrix<Ty> m3(dim, dim, '0');*/		// in release mode, int, GF2, GF2_test is same, thanks to the compiler

		start = clock();	//start timing
		for (int i = 0; i < 1; ++i)
			m3 = m1 * m2;		// with my_double, in debug mode 3.7s, in release mode, 2.9s, this is faster for POD

		end = clock();		//end timing
		cout.setf(ios::scientific);
		cout << scientific << "time consume: " << ((double)end - start) / CLOCKS_PER_SEC << "s\n\n";
		cout.unsetf(ios::floatfield);
	}
	static void matrix_inv() {

		const int n = 6;
		Matrix<Complex> A(n, n);
		int n2 = n * n;
		for (int i = 0; i < n2; ++i) {
			A(i) = Complex(my::rand_u(), my::rand_u());
		}
		Matrix<Complex> A_inv = A.inv();
		//cout << "A_inv" << A_inv << endl;
#ifdef count_operation_number
		cout << "operation_number=" << my_double_auxiliary_storage::operation_number << endl;
#endif
		cout << "A*A_inv" << A * A_inv;

#ifdef count_operation_number
		cout << "operation_number=" << my_double_auxiliary_storage::operation_number << endl;
#endif
	}
	static void matrix_end_max() {
		Matrix<int> A(6, 3, my::rand_int);
		cout << "A" << A;
		int pos;
		int me = A.max_abs_ele_end_from(3, pos);
		cout << "pos=" << pos << endl;
		cout << "me=" << me << endl;

		me = A.max_abs_col_ele_end_from(0, 1, pos);
		cout << "pos=" << pos << endl;
		cout << "me=" << me << endl;

		A = Matrix<int>(6, 3, 'i');
		cout << "A" << A;
		me = A.non_0_col_ele_end_from(0, 5, pos);
		cout << "pos=" << pos << endl;
		cout << "me=" << me << endl;

	}
	static void matrix_low_triangle() {
		const int m = 6;
		const int n = 8;
		Matrix<my_double> A(m, n, {
			12, 35, 1, 6,26,19,24,13,
			50, 3,32, 7,21,23,25,10.5,
			39, 31, 9, 2,22,27,20,11,
			61, 8,28,33,17,10,15,8.5,
			1, 30, 5,34,12,14,16,6,
			7.4, 4,36,29,13,18,11,6.5
			});
		cout << "A" << A;
		int p = A.row_transformation_to_low_triangle();
		cout << "A" << A;

		Matrix<int> permute_record = A.col_permute_to_full_rank_on_right();
		cout << "A" << A;
		cout << "permute_record" << permute_record;
		A.row_transformation_right_low_triangle_to_identity();
		cout << "A" << A;
	}
	static void matrix_up_triangle() {
		const int m = 6;
		const int n = 8;
		Matrix<my_double> A(m, n, {
			12, 35, 1, 6,26,19,24,13,
			50, 3,32, 7,21,23,25,10.5,
			39, 31, 9, 2,22,27,20,11,
			61, 8,28,33,17,10,15,8.5,
			1, 30, 5,34,12,14,16,6,
			7.4, 4,36,29,13,18,11,6.5
			});
		A.switch_col(0, 7);
		cout << "A" << A;
		A.row_transformation_to_up_triangle();
		cout << "A" << A;

		Matrix<int> permute_record = A.col_permute_to_full_rank_on_left();
		cout << "A" << A;
		cout << "permute_record" << permute_record;
		A.row_transformation_left_up_triangle_to_identity();
		cout << "A" << A;
	}
	static void matrix_permute() {
		Matrix<int> A(1, 6, { 1,0,3,5,2,4 });
		Matrix<int> B(3, 6, my::rand_int);
		cout << "A" << A;
		cout << "B" << B;
		B.permute_col(A);
		cout << "B" << B;
		Matrix<int> C(1, 6, { 10,20,30,40,50,60 });
		cout << "C" << C;
		C.permute(A);
		cout << "C" << C;
		C.permute_back(A);
		cout << "C" << C;
	}

	/* bch code */
	static void bch_parity_matrix() {

		const int m = 5;
		const int t = 2;
		GF2e<m>::init();

		BCH<m, t> bch;
		Matrix<GF2> pm = bch.get_parity_matrix();
		bch.print_info();
		cout << "bch.get_parity_matrix()" << pm;
		cout << "bch.get_parity_matrix() * bch.get_generator_matrix().Transpose()" << pm * bch.get_generator_matrix().Transpose();
	}
	static void bch_info() {

		const int m = 6;	// GF field of 2^m
		const int t = 4;	// error that can correct
		GF2e<m>::init();
		BCH<m, t> bch;
		bch.print_info();
	}
	static void bch_simulate() {

		const int m = 6;
		const int t = 4;		// t=2, (63,51,5); t=3, (63,45,7); t=4, (63,39,9)
		GF2e<m>::init();

		BCH<m, t> bch;
		decode_type dt = decode_type::OSD;
		const int eta = 2;
		const int order = 1;

		int iteration = 100;		// 3000000, 100
		const int SNR_num = 11;
		Matrix<my_double> SNR_db(1, SNR_num, { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 });

		bch.print_info();
		switch (dt) {
		case(decode_type::Algebra):cout << "Algebra simulation" << endl; break;
		case(decode_type::GMD):cout << "GMD simulation" << endl; break;
		case(decode_type::Chase2):cout << "Chase2 simulation, flip bit=" << eta << endl; break;
		case(decode_type::Chase3): cout << "Chase3 simulation" << endl; break;
		case(decode_type::OSD):cout << "OSD simulation, order=" << order << endl; break;
		case(decode_type::Hybrid_Chase2_OSD):cout << "Hybrid_Chase2_OSD simulation, order=" << order << ", eta=" << eta << endl; break;
		default: cout << "decode type input error" << endl;
		}

		cout << "error frames = " << 100 << ", Channel: AWGN, Modulation : BPSK" << endl;
		cout << "\n------------\n" << endl;

		Matrix<my_double> BER(1, SNR_num);
		Matrix<my_double> FER(1, SNR_num);
		Matrix<int> Frame(1, SNR_num);
		Matrix<my_double> GF2_operation(1, SNR_num);
		Matrix<my_double> GF2e_operation(1, SNR_num);
		Matrix<my_double> my_double_operation(1, SNR_num);		// simulation result

		clock_t start, end;
		start = clock();		// starting time

		for (int i = 0; i < SNR_num; i++) {

			simulation::linear_block_code(bch, SNR_db(i), iteration, simulation::_error_frame_upper_bound, dt, eta, order);
			BER(i) = simulation::BER;
			FER(i) = simulation::FER;
			Frame(i) = simulation::frame_proceed;
			GF2_operation(i) = simulation::GF2_operation_cnt_ave;
			GF2e_operation(i) = simulation::GF2e_operation_cnt_ave;
			my_double_operation(i) = simulation::my_double_operation_cnt_ave;
		}

		end = clock();

		cout.setf(ios::scientific);
		cout << scientific << "time consume: " << ((double)end - start) / CLOCKS_PER_SEC << "s\n\n";
		cout.unsetf(ios::floatfield);

		cout << "SNR_db" << SNR_db;

		cout.setf(ios::scientific);
		cout << "BER" << BER;
		cout.unsetf(ios::floatfield);

		cout.setf(ios::scientific);
		cout << "FER" << FER;
		cout.unsetf(ios::floatfield);

		cout << "Frame" << Frame;
		cout << "GF2_operation" << GF2_operation;
		cout << "GF2e_operation" << GF2e_operation;
		cout << "my_double_operation" << my_double_operation;
	}
	static void bch_simulate_auto() {

		const int m = 6;
		const int t = 4;		// t=2, (63,51,5); t=3, (63,45,7); t=4, (63,39,9)
		GF2e<m>::init();

		BCH<m, t> bch;
		decode_type dt = decode_type::Hybrid_Chase2_OSD;
		const int order = 1;
		const int eta = 4;

		int iteration = 100;		// 3000000, 100
		const int SNR_num = 30;			// large enougth SNR number
		Matrix<my_double> SNR_db(1, SNR_num, 'v');

		bch.print_info();
		switch (dt) {
		case(decode_type::Algebra):cout << "Algebra simulation" << endl; break;
		case(decode_type::GMD):cout << "GMD simulation" << endl; break;
		case(decode_type::Chase2):cout << "Chase2 simulation, flip bit=" << eta << endl; break;
		case(decode_type::Chase3): cout << "Chase3 simulation" << endl; break;
		case(decode_type::OSD):cout << "OSD simulation, order = " << order << endl; break;
		case(decode_type::Hybrid_Chase2_OSD):cout << "Hybrid_Chase2_OSD simulation, order = " << order << ", eta = " << eta << endl; break;
		default: cout << "decode type input error" << endl;
		}

		cout << "error frames = " << 100 << ", Channel: AWGN, Modulation : BPSK" << endl;
		cout << "\n------------\n" << endl;

		Matrix<my_double> BER(1, SNR_num, 'v');
		Matrix<my_double> FER(1, SNR_num, 'v');
		Matrix<int> Frame(1, SNR_num, 'v');
		Matrix<my_double> GF2_operation(1, SNR_num, 'v');
		Matrix<my_double> GF2e_operation(1, SNR_num, 'v');
		Matrix<my_double> my_double_operation(1, SNR_num, 'v');		// simulation result

		clock_t start, end;
		start = clock();		// starting time

		SNR_db.push_back(0.0);
		do {
			simulation::linear_block_code(bch, SNR_db.back(), iteration, simulation::_error_frame_upper_bound, dt, eta, order);
			BER.push_back(simulation::BER);
			FER.push_back(simulation::FER);
			Frame.push_back(simulation::frame_proceed);
			GF2_operation.push_back(simulation::GF2_operation_cnt_ave);
			GF2e_operation.push_back(simulation::GF2e_operation_cnt_ave);
			my_double_operation.push_back(simulation::my_double_operation_cnt_ave);

			SNR_db.push_back(SNR_db.back() + 0.5);

		} while (Frame.back() != iteration);
		SNR_db.pop_back();

		end = clock();

		cout.setf(ios::scientific);
		cout << scientific << "time consume: " << ((double)end - start) / CLOCKS_PER_SEC << "s\n\n";
		cout.unsetf(ios::floatfield);

		cout << "SNR_db" << SNR_db;

		cout.setf(ios::scientific);
		cout << "BER" << BER;
		cout.unsetf(ios::floatfield);

		cout.setf(ios::scientific);
		cout << "FER" << FER;
		cout.unsetf(ios::floatfield);

		cout << "Frame" << Frame;
		cout << "GF2_operation" << GF2_operation;
		cout << "GF2e_operation" << GF2e_operation;
		cout << "my_double_operation" << my_double_operation;
	}

	/* ebch code */
	static void ebch_parity_matrix() {

		const int m = 7;
		const int t = 10;
		GF2e<m>::init();

		eBCH<m, t> ebch;
		Matrix<GF2> pm = ebch.get_parity_matrix();
		ebch.print_info();
		cout << "ebch.get_parity_matrix()" << pm;
		cout << "ebch.get_parity_matrix() * ebch.get_generator_matrix().Transpose()" << pm * ebch.get_generator_matrix().Transpose();
	}
	static void ebch_info() {

		const int m = 4;	// GF field of 2^m
		const int t = 2;	// error that can correct
		GF2e<m>::init();
		eBCH<m, t> ebch;
		ebch.print_info();
	}

	/* simulation function*/
	static vector<simulation_param> simulation_param_set_general;
	static void run_operation() {
		GF2e<test_OSD_Chase_orig::_m_>::init();
		BCH<test_OSD_Chase_orig::_m_, test_OSD_Chase_orig::_t_> bch;
		simulate_operation(bch);
	}
	static void run_FER_auto() {
		GF2e<test_OSD_Chase_orig::_m_>::init();
		BCH<test_OSD_Chase_orig::_m_, test_OSD_Chase_orig::_t_> bch;
		simulate_FER_auto(bch, 0.0, 0.5);
	}
	static void run_FER_specified_SNR() {
		GF2e<test_OSD_Chase_orig::_m_>::init();
		BCH<test_OSD_Chase_orig::_m_, test_OSD_Chase_orig::_t_> bch;
		Matrix<my_double> SNR_db(1, 2, { 5,5.2 });
		simulate_FER_specified_SNR(bch, SNR_db);
	}
	static void run_operation_all_in_one() {
		/*------------------n=63-----------------------*/
		GF2e<test_OSD_Chase_orig::_m1_>::init();
		BCH<test_OSD_Chase_orig::_m1_, test_OSD_Chase_orig::_t11_> bch_11;
		simulation_param_set_general = { {decode_type::Chase2,10},{decode_type::Chase2,12} };
		simulate_operation(bch_11);

		BCH<test_OSD_Chase_orig::_m1_, test_OSD_Chase_orig::_t12_> bch_12;
		simulation_param_set_general = { {decode_type::Chase2,8},{decode_type::Chase2,10},{decode_type::Chase2,12} };
		simulate_operation(bch_12);

		BCH<test_OSD_Chase_orig::_m1_, test_OSD_Chase_orig::_t13_> bch_13;
		simulation_param_set_general = { {decode_type::Hybrid_Chase2_OSD,1,2},{decode_type::Hybrid_Chase2_OSD,1,4} };
		simulate_operation(bch_13);

		/*------------------n=127-----------------------*/
		GF2e<test_OSD_Chase_orig::_m2_>::init();
		BCH<test_OSD_Chase_orig::_m2_, test_OSD_Chase_orig::_t21_> bch_21;
		simulation_param_set_general = { {decode_type::Hybrid_Chase2_OSD,1,10},{decode_type::Hybrid_Chase2_OSD,1,12}, \
		{decode_type::Hybrid_Chase2_OSD, 2, 4}, { decode_type::Hybrid_Chase2_OSD,2,8 }, { decode_type::Chase2,10 }, \
			{ decode_type::Chase2,12 }, { decode_type::OSD,3 }, };
		simulate_operation(bch_21);

		/*------------------n=255-----------------------*/
		GF2e<test_OSD_Chase_orig::_m3_>::init();
		BCH<test_OSD_Chase_orig::_m3_, test_OSD_Chase_orig::_t31_> bch_31;
		simulation_param_set_general = { {decode_type::Hybrid_Chase2_OSD,1,10},{decode_type::Hybrid_Chase2_OSD,1,12}, \
		{decode_type::Hybrid_Chase2_OSD, 2, 8}, { decode_type::Hybrid_Chase2_OSD,2,12 }, { decode_type::Chase2,10 }, \
		{ decode_type::Chase2, 12 }, {decode_type::OSD,3}, };
		simulate_operation(bch_31);

		BCH<test_OSD_Chase_orig::_m3_, test_OSD_Chase_orig::_t32_> bch_32;
		simulation_param_set_general = { { decode_type::Chase2,10 }, { decode_type::Chase2,12 }, \
		{ decode_type::Hybrid_Chase2_OSD,1,10 }, { decode_type::Hybrid_Chase2_OSD,1,12 }, { decode_type::Hybrid_Chase2_OSD,2,10 }, };
		simulate_operation(bch_32);

		BCH<test_OSD_Chase_orig::_m3_, test_OSD_Chase_orig::_t33_> bch_33;
		simulation_param_set_general = { { decode_type::Chase2,10 }, { decode_type::Chase2,12 }, };
		simulate_operation(bch_33);
	}
	static void compute_OSD_dependent_col() {
		GF2e<test_OSD_Chase_orig::_m_>::init();
		eBCH<test_OSD_Chase_orig::_m_, test_OSD_Chase_orig::_t_> ebch;
		Matrix<my_double> SNR_db = (my_double)4.0;		// it doesn't matter what SNR takes
		//test_OSD_Chase_orig::simulate_FER_specified_SNR(bch, SNR_db);
		simulation::_iteration = 30000;
		simulation::_error_frame_upper_bound = 30000;

		simulate_param(ebch, { decode_type::OSD,1 }, SNR_db);

		cout << "OSD::sum_dependent_column" << OSD::sum_dependent_column;
	}
};

vector<simulation_param> test_OSD_Chase_orig::simulation_param_set_general = {

	//basic parameter
	/*{decode_type::OSD,2},
	{decode_type::OSD,1},
	{decode_type::Chase2,2},
	{decode_type::Chase2,4},
	{decode_type::Chase2,6},
	{decode_type::Chase2,8},
	{decode_type::Hybrid_Chase2_OSD,1,6},
	{decode_type::Hybrid_Chase2_OSD,1,8},
	{decode_type::Algebra},
	{decode_type::GMD},*/

	// additional parameter
	{decode_type::Chase2,7},
	{decode_type::Chase2,10},
};
