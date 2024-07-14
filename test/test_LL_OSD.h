/*****************************************************************//**
 * \file   LC_OSD_test.h
 * \brief  pick up LL OSD
 *
 * \author 26259
 * \date   March 2024
 *********************************************************************/

#pragma once

#include"test_common.h"

extern int frame_max;
extern int frame_min;
extern int frame_error_max;
extern my_double SNR_dB_start;
extern my_double SNR_dB_step;
extern my_double SNR_dB_end;

class test_LL_OSD {
public:

	static bool if_output_to_file;

	static void init_frame_SNR() {

		bool run_performance = true;
		if_output_to_file = false;

		if (run_performance) {
			frame_max = 50000000;		// 50000000
			frame_min = 10000;			// 10000
			frame_error_max = 200;		// 200
			SNR_dB_start = 3;		// 0
			SNR_dB_step = 1;		// 1,  0.5 
			SNR_dB_end = 8;		// 8, 7, 4.5, 4
		}
		else {
			frame_max = 100000;
			frame_min = 100000;
			frame_error_max = 1000000;
			SNR_dB_start = 4;		// 0
			SNR_dB_step = 1;		// 1,  0.5 
			SNR_dB_end = 6;		// 8, 7, 4.5, 4
		}
	}

	static int LL_OSD_FER_old() {		// incorrect simulation, higher FER than normal
		init_frame_SNR();
		//cout << boolalpha;

		// BCH code parameter
		const int m = 6;
		const int t = 3;
		// LL_OSD parameter
		const int order = 3;

		//// BCH code parameter
		//const int m = 6;
		//const int t = 3;
		//// LC_OSD parameter
		//const int delta = 6;
		//const int Lmax = 1024;

		GF2e<m>::init();
		BCH<m, t> bch;
		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();
		const int n = (1 << m) - 1;
		int k = bch.get_k();
		const int d = 2 * t + 1;
		const int k_prime = n - d + 1;

		/*cout << "G = " << endl;
		G.print();
		cout << "H = " << endl;
		H.print();*/

		cout << "-------------------------" << endl << endl;

		Matrix<GF2> u(1, k);											// message
		for (int i = 0; i < k; ++i) {
			u[i] = my::rand_01();
		}
		Matrix<GF2> c = u * G;									// transmitted codeword
		cout << "c = " << endl;
		c.print();
		cout << "-------------------------" << endl << endl;

		cout << "BCH(" << n << ", " << k << ", " << d << ")" << endl;


		Matrix<my_double> SNR_dB_vec(SNR_dB_start, SNR_dB_step, SNR_dB_end, 'd');
		//Vector_ext::equally_distant(SNR_dB_start, SNR_dB_end, SNR_dB_step, SNR_dB_vec);// 	SNR_dB_vec[0] = 6;
		//cout << "SNR_dB_vec" << SNR_dB_vec;


		LL_OSD<m, t, k_prime> ll_osd;
		cout << "LL_OSD(" << order << ")" << endl;

		// set cout into file
		if (if_output_to_file) {
			char file_name[155] = { 0 };
			sprintf_s(file_name, 155, "(%d,%d,%d)BCH_LL_OSD(%d).txt", n, k, d, order);
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}
		else {
			printf("\n");
		}
		printf("SNR(dB)   FER          Errors    Frames        GF2_ope      GF2e_ope     float_ope    ave_list     invalid_list_rate      ET_rate          time(s)\n");
		clock_t start, end;

		my::set_seed_adv(1);

		for (int i = 0; i < SNR_dB_vec.size(); ++i) {

			double SNR = (double)pow(10, SNR_dB_vec[i] / 10);
			double sigma = sqrt(n / (2 * k * SNR));
			int frame_error = 0;
			int frame_finished = 0;
			GF2_auxiliary_storage::operation_number = 0;
			my_double_auxiliary_storage::operation_number = 0; 
			GF2e_auxiliary_storage::operation_number = 0;

			ll_osd.total_used_list_num = 0;
			ll_osd.num_invalid_list = 0;
			ll_osd.num_early_stop = 0;

			start = clock();
			while (frame_finished < frame_max) {

				/*Matrix<my_double, 1, n> r({
					-0.929966,-0.777951,-0.852718,1.35274,0.0624572,-0.267671,-0.690881,1.07161,-0.37672,-1.59395,-0.84877,-1.48822,-1.39876,-0.137612,0.991339,
				});*/										// recieved word


				Matrix<my_double> r(1, n);
				for (int p = 0; p < n; ++p) {
					r[p] = (c[p] == 0 ? 1 : -1) + my::rand_ga_adv() * (double)sigma;	// BPSK modulation and AWGN channel
				}


				Matrix<GF2> c_hat = ll_osd.decode_v(r, order);

				if (c_hat != c) {
					frame_error++;
					//cout << "decode error r";
					//r.print();
				}

				frame_finished++;
				if (frame_error >= frame_error_max && frame_finished >= frame_min) {
					break;
				}
			}

			end = clock();

			double time_consume = ((double)end - start) / CLOCKS_PER_SEC;
			//cout << "time_consume = " << time_consume << " s" << endl;

			printf("%2.1f       %0.3E    %5d     %10d    %0.3E    %0.3E    %0.3E    %0.3E    %0.3E              %0.3E        %0.3E\n", \
				(double)SNR_dB_vec[i], (double)frame_error / frame_finished, frame_error, frame_finished, \
				(double)GF2_auxiliary_storage::operation_number / frame_finished, \
				(double)GF2e_auxiliary_storage::operation_number / frame_finished, \
				(double)my_double_auxiliary_storage::operation_number / frame_finished, \
				(double)ll_osd.total_used_list_num / frame_finished, \
				(double)ll_osd.num_invalid_list / frame_finished, \
				(double)ll_osd.num_early_stop / frame_finished, \
				time_consume / frame_finished);
		}

		return 0;
	}

	static int LL_OSD_FER() {

		init_frame_SNR();
		//cout << boolalpha;

		// BCH code parameter
		const int m = 5;
		const int t = 1;

		// LL_OSD parameter
		int order = 1;
		int second_order = 0;		// in segmentation mode, ensure order + second_order == max order is enough

		// LL_OSD Hybrid Viterbi parameter
		int order_hybrid = 1, selected_row_num = 6, max_list_num = 128;

		GF2e<m>::init();
		BCH<m, t> bch;
		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();
		int n = G.col();
		int k = G.row();
		int d = bch.get_d();

		/*cout << "G = " << endl;
		G.print();
		cout << "H = " << endl;
		H.print();*/

		cout << "-------------------------" << endl << endl;

		Matrix<GF2> u(1, k);											// message
		for (int i = 0; i < k; ++i) {
			u[i] = my::rand_01();
		}
		Matrix<GF2> c = u * G;									// transmitted codeword
		cout << "c = " << endl;
		c.print();
		cout << "-------------------------" << endl << endl;

		cout << "BCH(" << n << ", " << k << ", " << d << ")" << endl;

		Matrix<my_double> SNR_dB_vec(SNR_dB_start, SNR_dB_step, SNR_dB_end, 'd');
		//Vector_ext::equally_distant(SNR_dB_start, SNR_dB_end, SNR_dB_step, SNR_dB_vec);// 	SNR_dB_vec[0] = 6;
		//cout << "SNR_dB_vec" << SNR_dB_vec;

		LL_OSD_v2<m,t> ll_osd(order);		//osd.ibuc.reset_s(s);
		//LL_OSD_binary_H_v2<m, t> ll_osd(order, second_order);		//osd.ibuc.reset_s(s);
		
		//LL_OSD_binary_H_hybrid_Viterbi_v2<m, t> ll_osd(order_hybrid, selected_row_num, max_list_num);		//osd.ibuc.reset_s(s);

		/*ll_osd.seg_weight_sequence = Matrix<int>(6, 2, {
			0,1,
			1,0,
			0,2,
			1,1,
			0,3,
			1,2
		});
		cout << "ll_osd.seg_weight_sequence" << ll_osd.seg_weight_sequence;*/

		ll_osd.ET_type = 4;			// restrict to ML criterion in segmented LLOSD

		// set cout into file
		if (if_output_to_file) {
			char file_name[155] = { 0 };
			if (ll_osd.type == 10)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_LL_OSD_v2(%d).txt", n, k, d, order);		// it semms harder and no state decrease.
			else if (ll_osd.type == 11)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_LL_OSD_v2_seg(%d,%d).txt", n, k, d, order, second_order);		// it semms harder and no state decrease.
			else if (ll_osd.type == 12)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_LL_OSD_binary_H(%d).txt", n, k, d, order);		// it semms harder and no state decrease.
			else if (ll_osd.type == 13)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_LL_OSD_binary_H_seg(%d,%d).txt", n, k, d, order, second_order);		// it semms harder and no state decrease.
			else if (ll_osd.type == 14)
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_LL_OSD_binary_H_Hybrid_Viterbi(%d,%d,%d).txt", n, k, d, order_hybrid, selected_row_num, max_list_num);
			else
				sprintf_s(file_name, 155, "(%d,%d,%d)BCH_OSD_unknown(%d).txt", n, k, d, order);
			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}
		else {
			printf("\n");
		}
		//printf("SNR(dB)   FER          Errors    Frames        GF2_ope      GF2e_ope     float_ope    time(s)\n");
		printf("SNR(dB)   FER          Errors    Frames        GF2_ope      GF2e_ope     float_ope    ave_list     binary_list    ET_rate        invalid_rate      time(s)\n");
		//printf("________________________________________________________________________________\n\n");
		clock_t start, end;

		my::set_seed_adv(1);

		for (int i = 0; i < SNR_dB_vec.size(); ++i) {

			double SNR = (double)pow(10, SNR_dB_vec[i] / 10);
			double sigma = sqrt(n / (2 * k * SNR));
			int frame_error = 0;
			int frame_finished = 0;
			GF2_auxiliary_storage::operation_number = 0;
			GF2e_auxiliary_storage::operation_number = 0;
			my_double_auxiliary_storage::operation_number = 0;
			unsigned long long TEP_total = 0;
			unsigned long long binary_list_total = 0;
			int ET_frame_total = 0;
			int invalid_frame_total = 0;

			ll_osd.sigma = sigma;


			start = clock();
			while (frame_finished < frame_max) {

				//Matrix<my_double> r(1, n, {
				//	-0.4137,1.52826,-1.83037,0.937352,-0.0571138,0.3476,0.909213,
				//});										// recieved word

				Matrix<my_double> r(1, n);
				for (int p = 0; p < n; ++p) {
					r[p] = ((int)c[p] == 0 ? 1 : -1) + my::rand_ga_adv() * (double)sigma;	// BPSK modulation and AWGN channel
				}

				//Vector<int, n> nat;
				//Vector_ext::natual(nat);
				//cout << "nat" << nat;
				//cout << "r";
				//r.print();

				//ll_osd.call_GE_only(r);			// test GE complexity
				ll_osd.solve(r);				// estimated codeword
				if (ll_osd.c_hat != c) {
					frame_error++;

					/*cout << "frame_finished = " << frame_finished << endl;
					cout << "decode error r\n";
					r.print();*/

					/*cout << "ll_osd.c_hat\n"; 
					ll_osd.c_hat.print();*/
				}
				TEP_total += ll_osd.TEP_num;
				binary_list_total += ll_osd.binary_codeword_num;
				ET_frame_total += (ll_osd.is_early_termination == true);
				invalid_frame_total += (ll_osd.has_codeword_binary == false);

				//osd_dual.solve(r);				// estimated codeword
				//if (osd_dual.c_hat != c) {
				//	frame_error++;
				//	cout << "decode error r";
				//	r.print();
				//	osd.solve(r);				// estimated codeword
				//	if (osd.c_hat == c) {					
				//		cout << "decode should be correct" << endl;
				//	}
				//}

				frame_finished++;
				if (frame_error >= frame_error_max && frame_finished >= frame_min) {
					break;
				}
			}

			end = clock();

			double time_consume = ((double)end - start) / CLOCKS_PER_SEC;
			//cout << "time_consume = " << time_consume << " s" << endl;

			//printf("SNR_dB= %2.1f ;FER= %0.3E ;Errors= %5d ;Frames= %10d ;GF2_ope= %0.3E ;float_ope= %0.3E ;time(s)= %0.3E\n", \
			//	(double)SNR_dB_vec[i], (double)frame_error / frame_finished, frame_error, frame_finished, \
			//	(double)GF2_auxiliary_storage::operation_number / frame_finished, \
			//	(double)my_double_auxiliary_storage::operation_number / frame_finished, time_consume / frame_finished);
			printf("%2.1f       %0.3E    %5d     %10d    %0.3E    %0.3E    %0.3E    %0.3E    %0.3E      %0.3E      %0.3E         %0.3E\n", \
				(double)SNR_dB_vec[i], (double)frame_error / frame_finished, frame_error, frame_finished, \
				(double)GF2_auxiliary_storage::operation_number / frame_finished, \
				(double)GF2e_auxiliary_storage::operation_number / frame_finished, \
				(double)my_double_auxiliary_storage::operation_number / frame_finished, \
				(double)TEP_total / frame_finished, \
				(double)binary_list_total / frame_finished, \
				(double)ET_frame_total / frame_finished, \
				(double)invalid_frame_total / frame_finished, \
				time_consume / frame_finished);
		}

		return 0;
	}

	static int LL_OSD_Viterbi_FER() {

		init_frame_SNR();
		//cout << boolalpha;

		// BCH code parameter
		const int m = 3;
		const int t = 1;
		const int n = (1 << m) - 1;
		const int d = 2 * t + 1;
		const int k_prime = n - d + 1;

		/* decoding parameter */
		const int selected_row_num = 2;
		const int max_list_num = 4;

		GF2e<m>::init();
		BCH<m, t> bch;
		Matrix<GF2> G = bch.get_generator_matrix();
		Matrix<GF2> H = bch.get_parity_matrix();
		int k = G.row();

		/*cout << "G = " << endl;
		G.print();
		cout << "H = " << endl;
		H.print();*/

		cout << "-------------------------" << endl << endl;

		Matrix<GF2> u(1, k);											// message
		for (int i = 0; i < k; ++i) {
			u[i] = my::rand_01();
		}
		Matrix<GF2> c = u * G;									// transmitted codeword
		cout << "c = " << endl;
		c.print();
		cout << "-------------------------" << endl << endl;

		cout << "BCH(" << n << ", " << k << ", " << d << ")" << endl;

		Matrix<my_double> SNR_dB_vec(SNR_dB_start, SNR_dB_step, SNR_dB_end, 'd');
		//Vector_ext::equally_distant(SNR_dB_start, SNR_dB_end, SNR_dB_step, SNR_dB_vec);// 	SNR_dB_vec[0] = 6;
		//cout << "SNR_dB_vec" << SNR_dB_vec;

		LL_OSD_Viterbi<m, t, k_prime, selected_row_num> LL_vit;

		// set cout into file
		if (if_output_to_file) {
			char file_name[155] = { 0 };
			sprintf_s(file_name, 155, "(%d,%d,%d)BCH_LL_OSD_Viterbi(%d,%d).txt", n, k, d, selected_row_num, max_list_num);		// it semms harder and no state decrease.

			FILE* stream1;
			freopen_s(&stream1, file_name, "w", stdout);
		}
		else {
			printf("\n");
		}
		printf("SNR(dB)   FER          Errors    Frames        GF2_ope      GF2e_ope     float_ope    ave_list     ET_rate        invalid_rate      time(s)\n");
		//printf("________________________________________________________________________________\n\n");
		clock_t start, end;

		my::set_seed_adv(1);

		for (int i = 0; i < SNR_dB_vec.size(); ++i) {

			double SNR = (double)pow(10, SNR_dB_vec[i] / 10);
			double sigma = sqrt(n / (2 * k * SNR));
			int frame_error = 0;
			int frame_finished = 0;
			GF2_auxiliary_storage::operation_number = 0;
			GF2e_auxiliary_storage::operation_number = 0;
			my_double_auxiliary_storage::operation_number = 0;

			LL_vit.total_used_list_num = 0;
			LL_vit.num_invalid_list = 0;
			LL_vit.num_early_stop = 0;

			start = clock();
			while (frame_finished < frame_max) {


				//Matrix<my_double> r(1, n, {
				//	-0.4137,1.52826,-1.83037,0.937352,-0.0571138,0.3476,0.909213,
				//	});										// recieved word


				Matrix<my_double> r(1, n);
				for (int p = 0; p < n; ++p) {
					r[p] = (c[p] == 0 ? 1 : -1) + my::rand_ga_adv() * (double)sigma;	// BPSK modulation and AWGN channel
				}


				//Vector<int, n> nat;
				//Vector_ext::natual(nat);
				//cout << "nat" << nat;
				//cout << "r";
				//r.print();


				Matrix<GF2> c_hat = LL_vit.decode_v(r, max_list_num);				// estimated codeword
				if (c_hat != c) {
					frame_error++;

					/*cout << "frame_finished = " << frame_finished << endl;
					cout << "decode error r\n";
					r.print();*/

				}

				//osd_dual.solve(r);				// estimated codeword
				//if (osd_dual.c_hat != c) {
				//	frame_error++;
				//	cout << "decode error r";
				//	r.print();
				//	osd.solve(r);				// estimated codeword
				//	if (osd.c_hat == c) {					
				//		cout << "decode should be correct" << endl;
				//	}
				//}

				frame_finished++;
				if (frame_error >= frame_error_max && frame_finished >= frame_min) {
					break;
				}
			}

			end = clock();

			double time_consume = ((double)end - start) / CLOCKS_PER_SEC;
			//cout << "time_consume = " << time_consume << " s" << endl;

			//printf("SNR_dB= %2.1f ;FER= %0.3E ;Errors= %5d ;Frames= %10d ;GF2_ope= %0.3E ;float_ope= %0.3E ;time(s)= %0.3E\n", \
			//	(double)SNR_dB_vec[i], (double)frame_error / frame_finished, frame_error, frame_finished, \
			//	(double)GF2_auxiliary_storage::operation_number / frame_finished, \
			//	(double)my_double_auxiliary_storage::operation_number / frame_finished, time_consume / frame_finished);

			printf("%2.1f       %0.3E    %5d     %10d    %0.3E    %0.3E    %0.3E    %0.3E    %0.3E      %0.3E         %0.3E\n", \
				(double)SNR_dB_vec[i], (double)frame_error / frame_finished, frame_error, frame_finished, \
				(double)GF2_auxiliary_storage::operation_number / frame_finished, \
				(double)GF2e_auxiliary_storage::operation_number / frame_finished, \
				(double)my_double_auxiliary_storage::operation_number / frame_finished, \
				(double)LL_vit.total_used_list_num / frame_finished, \
				(double)LL_vit.num_early_stop / frame_finished, \
				(double)LL_vit.num_invalid_list / frame_finished, \
				time_consume / frame_finished);
		}

		return 0;
	}

};

bool test_LL_OSD::if_output_to_file = false;