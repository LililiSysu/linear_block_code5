/*****************************************************************//**
 * \file   Polar.h
 * \brief  adding the method of polar codes, developed by April 2022, this is still unused
 * 
 * \author lilili
 * \date   April 2023
 *********************************************************************/

// this file is independent to others in the project

#pragma once
#include<stdio.h>
#include<vector>
#include<algorithm>
#include<random>
#include<chrono>
#include<cmath>
#include<ctime>
using namespace std;
using namespace chrono;

class polar {
	static const unsigned short crc_table[256];

private:
	// private member

	int n;						// number of layer in polar code
	int N;						// number of bit (or channel) in polar code

	vector<int> u;				// the vector of information bit to encode, with size less than 2^n

	vector<int> v;				// the vector of encoding result, with size equals 2^n
	vector<int> v_sys;			// the vector of systematic encoding result, with size equals 2^n

	bool SNR_is_set;			// whether 'SNR' has set
	double SNR;					// Signal Noise Ratio in dB, that is Eb / N0 where Eb = 1, N0 = 2 * sigma^2
	double sigma;				// sigma value of Gaussian noise
	vector<double> received_v;	// 'v' after BPSK modulation in AWGN channel
	vector<int> decoded_u;		// the not freezed value of 'decoded_eu' after sorting by 'sorted_channel_ind'

	bool list_num_is_set;		// whether 'list_num' has set
	int list_num;				// number of path for list decoding

	bool rate_is_set;			// whether 'rate' is set
	double rate;				// rate of polar code, that is size of 'u' over 'N'

	int total_bit;				// number of information bit for error rate compute
	int total_SC_error_bit;		// number of error information bit for SC decode
	int total_fast_SC_error_bit;// number of error information bit for fast SC decode
	int total_SCL2_error_bit;	// number of error information bit for SCL deocde with list_num 2
	int total_SCL4_error_bit;	// number of error information bit for SCL deocde with list_num 4
	int total_SCL8_error_bit;	// number of error information bit for SCL deocde with list_num 8
	int total_SCL16_error_bit;	// number of error information bit for SCL deocde with list_num 16
	int total_SCL32_error_bit;	// number of error information bit for SCL deocde with list_num 32
	int total_SCL_error_bit;	// number of error information bit for SCL deocde with list_num not the above
	int total_fast_SCL_error_bit;// number of error information bit for SCL deocde with list_num not the above

	bool is_CRC16_aided;		// whether use CRC16 for polar code
	bool is_SCL_decode_failed;	// whether SCL decode failed
	int SCL_decode_failed_time;	// times of SCL decode failed
	int SCL_decode_failed_bit;	// number of bit cast due to SCL decode failed	
	bool use_adaptive_L;		// whether change L adaptively, see 'adaptive_fast_SCL_decode()'


	// private auxiliary member

	vector<pair<int, double>> sorted_channel_ind_ZWi;	// the vector to hold pair of index of ordered channel and its ZW value

	vector<double> LLR;									// the Log Likelilhood Ratio for each channel
	vector<int> active_channel_sum;						// the i^th element is the active channel number between index 0 and i-1
	vector<int> decoded_v_tmp;							// the v decoded from LLR (fast SC decode)
	vector<int> decoded_eu;								// the decoded result of 'received_v'

	vector<pair<int, double>> decoded_list_PM;			// the list of pair list index and decoded Path Metric, activate till size of list_num
	vector<vector<double>> decoded_list_LLR;			// list of decoded LLR, row equals 'list_num' should be configured first
	vector<vector<int>> decoded_list_eu;				// same as 'decoded_eu', but information bit stored in decoded_list_path
	vector<bool> occupied;								// whether the index of list occupied
	vector<bool> is_done;								// another indicator for list occupied

	double time_key_routine;							// time consume for key routine

	// private member function

	/**
	 * .default constructor function
	 *
	 */
	polar();

	/**
	 * .generate the vector to hold pair of index of ordered channel and its ZW value, changing 'sorted_channel_ind_ZWi'
	 *
	 */
	void generate_BEC_ZWi();

	/**
	 * .reverse bit order
	 *
	 * \param a:	number to solve
	 * \return		bit reverse of a
	 */
	int bit_revorder(int a);

	/**
	  * .polar encoding in recursive way, changing 'v'
	  *
	  * \param start_ind:	the start index of vector 'v' being processed in this recursion
	  * \param end_ind:		the end index of vector 'v' being processed in this recursion
	 */
	void encode_recur(int start_ind, int end_ind);

	/**
	  * .polar encoding in recursive way, changing 'v_sys'
	  *
	  * \param start_ind:	the start index of vector 'v_sys' being processed in this recursion
	  * \param end_ind:		the end index of vector 'v_sys' being processed in this recursion
	 */
	void encode_systematic_recur(int start_ind, int end_ind);

	/**
	 * .compute the upper part LLR in the former layer of polar code
	 *
	 * \param L1:	the first LLR
	 * \param L2:	the second LLR
	 * \return		the result LLR
	 */
	inline double f(double L1, double L2);

	/**
	 * .compute the lower part LLR in the former layer of polar code
	 *
	 * \param u:	the upper part encoded bit
	 * \param L1:	the upper in this layer LLR
	 * \param L2:	the lower in this layer LLR
	 * \return		the result LLR
	 */
	inline double g(int u, double L1, double L2);

	/**
	 * .get the minimum absolute value of x and y
	 *
	 * \param x:	the first input value
	 * \param y:	the second input value
	 * \return		result as stated
	 */
	inline double min_abs(double x, double y);

	/**
	 * .a simple two channel polar decode function, input as 'received_v', output as 'decode_u'
	 *
	 */
	void two_channel_decode();

	/**
	 * .print 'decoded_eu', cannot access now because in private member
	 *
	 */
	void print_decoded_eu();

	/**
	 * .decode polar code recursively, input as 'active_channel_sum', 'LLR', output as 'decoded_eu'
	 * changing 'LLR' and 'decoded_eu'
	 *
	 * \param LLR_start_ind:	start index of the part of LLR to decode
	 * \param LLR_end_ind:		end index of the part of LLR to decode
	 */
	void SC_decode_recur(int LLR_start_ind, int LLR_end_ind);

	/**
	 * .decode polar code recursively in a fast method, input as 'active_channel_sum', 'LLR',
	 * output as 'decoded_eu', changing 'LLR' and 'decoded_eu'
	 *
	 * \param LLR_start_ind:	start index of the part of LLR to decode
	 * \param LLR_end_ind:		end index of the part of LLR to decode
	 */
	void fast_SC_decode_recur(int LLR_start_ind, int LLR_end_ind);

	/**
	 * .polar encoding in recursive way, changing 'decoded_v_tmp'
	 *
	 * \param start_ind:	the start index of vector 'decoded_v_tmp' being processed in this recursion
	 * \param end_ind:		the end index of vector 'decoded_v_tmp' being processed in this recursion
	 */
	void encode_recur_for_dv(int start_ind, int end_ind);

	/**
	 * .compute the absolute value
	 *
	 * \param x:	number to solve
	 * \return		absolute value of 'x'
	 */
	inline double abs(double x);

	/**
	 * .list decode polar code recursively, input as 'active_channel_sum', 'decoded_list_LLR',
	 * 'decoded_list_eu', 'decoded_list_PM', output as 'decoded_eu'. not using lazy copy !!
	 *
	 * \param LLR_start_ind:	start index of the part of LLR to decode
	 * \param LLR_end_ind:		end index of the part of LLR to decode
	 */
	void SCL_decode_recur(int LLR_start_ind, int LLR_end_ind);

	/**
	 * .use LLR to store fixed received bit LLR to save copying bit
	 *
	 */
	void SCL_decode_recur_father();

	/**
	 * .fast list decode polar code recursively, only deal with R0 and REP case
	 *  input as 'active_channel_sum', 'decoded_list_LLR','decoded_list_eu',
	 * 'decoded_list_PM', output as 'decoded_eu'. not using lazy copy !!
	 *
	 * \param LLR_start_ind
	 * \param LLR_end_ind
	 */
	void fast_SCL_decode_recur(int LLR_start_ind, int LLR_end_ind);

	/**
	 * .use LLR to store fixed received bit LLR to save copying bit
	 *
	 */
	void fast_SCL_decode_recur_father();

public:
	// public constructor and destructor

	/**
	 * .construct with number of layer in polar code
	 *
	 * \param _n:	number of layer in polar code, greater than 0
	 */
	polar(int _n);

	/**
	 * .construct with number of layer in polar code and the vector of information bit to encode, with size less than 2^n
	 *
	 * \param _n:	number of layer in polar code, greater than 0
	 * \param _u:	the vector of information bit to encode, with size less than 2^n
	 */
	polar(int _n, const vector<int>& _u);

	// public member function

	/**
	 * .polar encode function, non-systematically, 'u' as input, changing 'v'
	 *
	 */
	void encode();

	/**
	 * .polar encode function, systematically, 'u' as input, changing 'v_sys'
	 *
	 */
	void encode_systematic();

	/**
	 * .change the vector of information bit to encode, 'v' and 'v_sys' is cleared
	 *
	 * \param new_u:	the new vector, with size less than 2^n,
	 */
	void set_u(const vector<int>& new_u);

	/**
	 * .print vector 'u'
	 *
	 */
	void print_u();

	/**
	 * .print vector 'v'
	 *
	 */
	void print_v();

	/**
	 * .print vector 'v_sys'
	 *
	 */
	void print_v_sys();

	/**
	 * .changing 'SNR' and 'sigma'
	 *
	 * \param _SNR:	the value to set SNR in dB
	 */
	void set_SNR(double _SNR);

	/**
	 * .input as 'v', changing 'received_v' using BPSK modulation, 0 -> 1 and 1 -> -1
	 *
	 */
	void modulation();

	/**
	 * .according to parameter 'sigma' to add Gaussian noise to 'received_v',
	 * should after modulation and set SNR!!
	 *
	 */
	void add_noise();

	/**
	 * .let the polar code pass the channel, input as 'v', output as 'received_v', should set SNR first !!
	 *
	 */
	void pass_channel();

	/**
	 * .SC decoding of polar code, input as 'received_v', 'sigma' and output as 'decoded_u'
	 *
	 */
	void SC_decode();

	/**
	 * .print 'decoded_u'
	 *
	 */
	void print_decoded_u();

	/**
	 * .fast SC decoding of polar code, input as 'received_v', 'sigma' and output as 'decoded_u'
	 *
	 */
	void fast_SC_decode();

	/**
	 * .set the list number for list decode
	 *
	 * \param _list_num:	the list number
	 */
	void set_list_num(int _list_num);

	/**
	 * .SC List decoding of polar code, input as 'received_v', 'sigma' and output as 'decoded_u'
	 *
	 */
	void SCL_decode();

	/**
	 * .set 'rate' for polar code
	 *
	 * \param _rate:	number to set
	 */
	void set_rate(double _rate);

	/**
	 * .randomly generete 'u' in {0,1}, 'rate' should be set first
	 *
	 */
	void generate_random_u();

	/**
	 * .count how many bit transmitted, result stored in 'total_bit'
	 *
	 */
	void information_bit_count();

	/**
	 * .count number of SC decode error bit, result stored in 'total_SC_error_bit'
	 *
	 */
	void SC_error_count();

	/**
	 * .count number of fast SC decode error bit, result stored in 'total_fast_SC_error_bit'
	 *
	 */
	void fast_SC_error_count();

	/**
	 * .count number of SCL decode error bit, result stored in 'total_SCL_error_bit'
	 *
	 */
	void SCL_error_count();

	/**
	 * .
	 *
	 * \return SC error rate, that is 'total_SC_error_bit' over 'total_bit'
	 */
	double SC_error_rate();

	/**
	 * .
	 *
	 * \return fast SC error rate, that is 'fast_total_SC_error_bit' over 'total_bit'
	 */
	double fast_SC_error_rate();

	/**
	 * .
	 *
	 * \return SCL error rate, that is 'total_SCL_error_bit' over 'total_bit'
	 */
	double SCL_error_rate();

	/**
	 * .clean error bit count for the next simulation
	 *
	 */
	void clear_error_count();

	/**
	 * .
	 *
	 * \param _is_CRC16_aided: true means set list decoding CRC16 aided
	 */
	void set_is_CRC16_aided(bool _is_CRC16_aided);

	/**
	 * .
	 *
	 * \return crc16 bit for 'u'
	 */
	unsigned short u_crc16();

	/**
	 * .put crc16 at back of 'u'
	 *
	 */
	void u_append_crc16();

	/**
	 * .remove the last 16 bit for 'u'
	 *
	 */
	void u_remove_crc16();

	/**
	 * .
	 *
	 * \return the last 16 bit in 'decoded_u'
	 */
	unsigned short du_get_crc16();

	/**
	 * .remove the last 16 bit for 'decoded_u'
	 *
	 */
	void du_remove_crc16();

	/**
	 * .
	 *
	 * \return crc16 value for 'decoded_u'
	 */
	unsigned short du_crc16();

	/**
	 * .append 16 bit at the end of 'decoded_u'
	 *
	 * \param us: bits to append
	 */
	void du_append_us(unsigned short us);

	/**
	 * .
	 *
	 * \return whether 'decoded_u' satisfied condition of CRC16
	 */
	bool is_du_satisfied_crc16();

	/**
	 * .
	 *
	 * \return 'time_key_routine'
	 */
	double get_time_key_routine();

	/**
	 * .fast SC List decoding of polar code, input as 'received_v', 'sigma' and output as 'decoded_u'
	 *
	 */
	void fast_SCL_decode();

	/**
	 * .count number of fast SCL decode error bit, result stored in 'total_fast_SCL_error_bit'
	 *
	 */
	void fast_SCL_error_count();

	/**
	 * .
	 *
	 * \return SCL error rate, that is 'total_fast_SCL_error_bit' over 'total_bit'
	 */
	double fast_SCL_error_rate();

	/**
	 * .adaptively increase parameter L as 1,2,4,8,16,32 until
	 * decode correction, if not at L 32, let the result be
	 * with minimum path metric
	 *
	 */
	void adaptive_fast_SCL_decode();
};


const unsigned short polar::crc_table[256] = {
0x0000, 0x1021, 0x2042, 0x3063, 0x4084, 0x50a5, 0x60c6, 0x70e7, 0x8108, \
0x9129, 0xa14a, 0xb16b, 0xc18c, 0xd1ad, 0xe1ce, 0xf1ef, 0x1231, 0x0210, \
0x3273, 0x2252, 0x52b5, 0x4294, 0x72f7, 0x62d6, 0x9339, 0x8318, 0xb37b, \
0xa35a, 0xd3bd, 0xc39c, 0xf3ff, 0xe3de, 0x2462, 0x3443, 0x0420, 0x1401, \
0x64e6, 0x74c7, 0x44a4, 0x5485, 0xa56a, 0xb54b, 0x8528, 0x9509, 0xe5ee, \
0xf5cf, 0xc5ac, 0xd58d, 0x3653, 0x2672, 0x1611, 0x0630, 0x76d7, 0x66f6, \
0x5695, 0x46b4, 0xb75b, 0xa77a, 0x9719, 0x8738, 0xf7df, 0xe7fe, 0xd79d, \
0xc7bc, 0x48c4, 0x58e5, 0x6886, 0x78a7, 0x0840, 0x1861, 0x2802, 0x3823, \
0xc9cc, 0xd9ed, 0xe98e, 0xf9af, 0x8948, 0x9969, 0xa90a, 0xb92b, 0x5af5, \
0x4ad4, 0x7ab7, 0x6a96, 0x1a71, 0x0a50, 0x3a33, 0x2a12, 0xdbfd, 0xcbdc, \
0xfbbf, 0xeb9e, 0x9b79, 0x8b58, 0xbb3b, 0xab1a, 0x6ca6, 0x7c87, 0x4ce4, \
0x5cc5, 0x2c22, 0x3c03, 0x0c60, 0x1c41, 0xedae, 0xfd8f, 0xcdec, 0xddcd, \
0xad2a, 0xbd0b, 0x8d68, 0x9d49, 0x7e97, 0x6eb6, 0x5ed5, 0x4ef4, 0x3e13, \
0x2e32, 0x1e51, 0x0e70, 0xff9f, 0xefbe, 0xdfdd, 0xcffc, 0xbf1b, 0xaf3a, \
0x9f59, 0x8f78, 0x9188, 0x81a9, 0xb1ca, 0xa1eb, 0xd10c, 0xc12d, 0xf14e, \
0xe16f, 0x1080, 0x00a1, 0x30c2, 0x20e3, 0x5004, 0x4025, 0x7046, 0x6067, \
0x83b9, 0x9398, 0xa3fb, 0xb3da, 0xc33d, 0xd31c, 0xe37f, 0xf35e, 0x02b1, \
0x1290, 0x22f3, 0x32d2, 0x4235, 0x5214, 0x6277, 0x7256, 0xb5ea, 0xa5cb, \
0x95a8, 0x8589, 0xf56e, 0xe54f, 0xd52c, 0xc50d, 0x34e2, 0x24c3, 0x14a0, \
0x0481, 0x7466, 0x6447, 0x5424, 0x4405, 0xa7db, 0xb7fa, 0x8799, 0x97b8, \
0xe75f, 0xf77e, 0xc71d, 0xd73c, 0x26d3, 0x36f2, 0x0691, 0x16b0, 0x6657, \
0x7676, 0x4615, 0x5634, 0xd94c, 0xc96d, 0xf90e, 0xe92f, 0x99c8, 0x89e9, \
0xb98a, 0xa9ab, 0x5844, 0x4865, 0x7806, 0x6827, 0x18c0, 0x08e1, 0x3882, \
0x28a3, 0xcb7d, 0xdb5c, 0xeb3f, 0xfb1e, 0x8bf9, 0x9bd8, 0xabbb, 0xbb9a, \
0x4a75, 0x5a54, 0x6a37, 0x7a16, 0x0af1, 0x1ad0, 0x2ab3, 0x3a92, 0xfd2e, \
0xed0f, 0xdd6c, 0xcd4d, 0xbdaa, 0xad8b, 0x9de8, 0x8dc9, 0x7c26, 0x6c07, \
0x5c64, 0x4c45, 0x3ca2, 0x2c83, 0x1ce0, 0x0cc1, 0xef1f, 0xff3e, 0xcf5d, \
0xdf7c, 0xaf9b, 0xbfba, 0x8fd9, 0x9ff8, 0x6e17, 0x7e36, 0x4e55, 0x5e74, \
0x2e93, 0x3eb2, 0x0ed1, 0x1ef0 };

// public constructor and destructor

polar::polar(int _n) :n(_n), N(1 << n) {
	generate_BEC_ZWi();
	SNR_is_set = false;
	list_num_is_set = false;
	rate_is_set = false;
	total_bit = 0;
	total_SC_error_bit = 0;
	total_fast_SC_error_bit = 0;
	total_SCL_error_bit = 0;
	is_CRC16_aided = false;
	is_SCL_decode_failed = false;
	SCL_decode_failed_time = 0;
	SCL_decode_failed_bit = 0;
	time_key_routine = 0;
	total_fast_SCL_error_bit = 0;
	use_adaptive_L = false;
}

polar::polar(int _n, const vector<int>& _u) : n(_n), N(1 << n), u(_u) {
	unsigned int ulen = u.size();
	if (ulen > (unsigned int)N)	throw "invalid 'u' size which greater than 2^n";
	generate_BEC_ZWi();
	active_channel_sum.resize(N + 1, 0);
	for (unsigned int i = 0; i < ulen; ++i) {
		active_channel_sum[sorted_channel_ind_ZWi[i].first + 1] = 1;
	}
	for (int i = 1; i <= N; ++i) {
		active_channel_sum[i] += active_channel_sum[i - 1];
	}
	SNR_is_set = false;
	list_num_is_set = false;
	rate_is_set = false;
	total_bit = 0;
	total_SC_error_bit = 0;
	total_fast_SC_error_bit = 0;
	total_SCL_error_bit = 0;
	is_CRC16_aided = false;
	is_SCL_decode_failed = false;
	SCL_decode_failed_time = 0;
	SCL_decode_failed_bit = 0;
	time_key_routine = 0;
	total_fast_SCL_error_bit = 0;
	use_adaptive_L = false;
}

// private member function

void polar::generate_BEC_ZWi() {
	sorted_channel_ind_ZWi.resize(N);
	for (int i = 0; i < N; ++i) {
		sorted_channel_ind_ZWi[i] = make_pair(i, 0);
	}
	sorted_channel_ind_ZWi[0].second = 0.5;	// initialize ZW to 0.5
	int m = 1;
	while (m <= N / 2) {
		for (int k = 0; k < m; ++k) {
			double z_tmp = sorted_channel_ind_ZWi[k].second;
			sorted_channel_ind_ZWi[k].second = 2 * z_tmp - z_tmp * z_tmp;	// use upper bound
			sorted_channel_ind_ZWi[k + m].second = z_tmp * z_tmp;
		}
		m <<= 1;
	}
	for (int i = 0; i < N; ++i) {
		sorted_channel_ind_ZWi[i].first = bit_revorder(sorted_channel_ind_ZWi[i].first);
	}
	sort(sorted_channel_ind_ZWi.begin(), sorted_channel_ind_ZWi.end(), \
		[](const pair<int, double>& a, const pair<int, double>& b) {
			return a.second < b.second;
		});
}

int polar::bit_revorder(int a) {
	int result = 0;
	for (int i = 0; a != 0; ++i) {
		if ((a & 1) == 0) {}
		else {
			result ^= 1 << (n - i - 1);
		}
		a >>= 1;
	}
	return result;
}

void polar::encode_recur(int start_ind, int end_ind) {
	if (start_ind == end_ind)	return;
	else {
		int len = end_ind - start_ind + 1;
		for (int i = start_ind; i < start_ind + len / 2; ++i) {
			v[i] = v[i] ^ v[i + len / 2];
		}
		encode_recur(start_ind, start_ind + len / 2 - 1);
		encode_recur(start_ind + len / 2, end_ind);
	}
}

void polar::encode() {
	v.resize(N);
	if (n == 0) {
		throw "error, case 'n'=0 is not a polar code";
	}
	for (int i = 0; i < N; ++i) {
		v[i] = 0;
	}
	int channel_num = u.size();

	for (int i = 0; i < channel_num; ++i) {
		v[sorted_channel_ind_ZWi[i].first] = u[i];	// use method of ZWi
	}
	encode_recur(0, v.size() - 1);
}

void polar::encode_systematic_recur(int start_ind, int end_ind) {
	if (start_ind == end_ind)	return;
	else {
		int len = end_ind - start_ind + 1;
		for (int i = start_ind; i < start_ind + len / 2; ++i) {
			v_sys[i] = v_sys[i] ^ v_sys[i + len / 2];
		}
		encode_systematic_recur(start_ind, start_ind + len / 2 - 1);
		encode_systematic_recur(start_ind + len / 2, end_ind);
	}
}

void polar::encode_systematic() {
	v_sys.resize(N);
	if (n == 0) {
		throw "error, case 'n'=0 is not a polar code";
	}
	for (int i = 0; i < N; ++i) {
		v_sys[i] = 0;
	}
	int channel_num = u.size();

	for (int i = 0; i < channel_num; ++i) {
		v_sys[sorted_channel_ind_ZWi[i].first] = u[i];			// use method of ZWi
	}
	encode_systematic_recur(0, v_sys.size() - 1);
	for (int i = N - 1; i >= channel_num; --i) {
		v_sys[sorted_channel_ind_ZWi[i].first] = 0;				// xAc=0, use method of ZWi
	}
	encode_systematic_recur(0, v_sys.size() - 1);
}

void polar::set_u(const vector<int>& new_u) {
	unsigned int ulen = u.size();
	u = new_u;
	if (is_CRC16_aided == true)	u_append_crc16();
	if ((int)u.size() > N)	throw "error of 'u' size greater than 'N'";
	if (ulen == u.size()) {}
	else {
		ulen = u.size();
		active_channel_sum.resize(N + 1, 0);
		for (unsigned int i = 0; i < ulen; ++i) {
			active_channel_sum[sorted_channel_ind_ZWi[i].first + 1] = 1;
		}
		for (int i = 1; i <= N; ++i) {
			active_channel_sum[i] += active_channel_sum[i - 1];
		}
	}
}

void polar::print_u() {
	int ulen;
	if (is_CRC16_aided == true)
		ulen = u.size() - 16;
	else
		ulen = u.size();
	printf("u    : ");
	for (int i = 0; i < ulen; ++i) {
		printf("%d ", u[i]);
	}
	printf("\n");
}

void polar::print_v() {
	if (v.size() != N)	throw "error size of 'v' not equals 'N'";
	printf("v    : ");
	for (int i = 0; i < N; ++i) {
		printf("%d ", v[i]);
	}
	printf("\n");
}

void polar::print_v_sys() {
	if (v_sys.size() != N)	throw "error size of 'v_sys' not equals 'N'";
	printf("v_sys: ");
	for (int i = 0; i < N; ++i) {
		printf("%d ", v_sys[i]);
	}
	printf("\n");
}

inline double polar::min_abs(double x, double y) {
	x = x < 0 ? -x : x;
	y = y < 0 ? -y : y;
	return x < y ? x : y;
}

inline double polar::f(double L1, double L2) {
	double m = min_abs(L1, L2);
	return ((L1 > 0 && L2 > 0) || (L1 < 0 && L2 < 0)) ? m : -m;
}

inline double polar::g(int u, double L1, double L2) {
	return (u == 0 ? L1 : -L1) + L2;
}

void polar::set_SNR(double _SNR) {
	if (rate_is_set == false)	throw "error of not setting rate";
	SNR = _SNR;
	//sigma = sqrt(0.5 * pow(10, -SNR / 10));
	sigma = sqrt(0.5 / rate * pow(10, -SNR / 10));	// use the same definition as book
	SNR_is_set = true;
}

void polar::modulation() {
	received_v.resize(N);
	for (int i = 0; i < N; ++i) {
		received_v[i] = v[i] == 0 ? 1 : -1;
	}
}

void polar::add_noise() {
	if (received_v.size() != N)	throw "error of 'received_v' size not equals to N";
	if (SNR_is_set == false)	throw "error of SNR not set, use function set_SNR to set SNR";
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<double> normal(0.0, sigma);

	for (int i = 0; i < N; ++i) {
		received_v[i] = received_v[i] + normal(gen);
	}
}

void polar::pass_channel() {
	modulation();
	add_noise();
}

void polar::two_channel_decode() {

	// to do: extend this to N channel, now done as function SC_decode

	decoded_eu.clear();
	// in this simple case vlen == 2
	LLR.resize(N);

	double scale = 2 / sigma / sigma;
	for (int i = 0; i < N; ++i) {
		LLR[i] = received_v[i] * scale;
	}

	// channel 0, update LLR
	LLR.push_back(f(LLR[0], LLR[1]));
	if (active_channel_sum[0 + 1] != active_channel_sum[0]) {
		// use LLR value
		decoded_eu.push_back(LLR.back() < 0);
	}
	else {
		// use freezed value
		decoded_eu.push_back(0);
	}
	LLR.pop_back();

	// channel 1, update LLR
	LLR.push_back(g(decoded_eu[0], LLR[LLR.size() - 2], LLR[LLR.size() - 1]));
	if (active_channel_sum[1 + 1] != active_channel_sum[1]) {
		// use LLR value
		decoded_eu.push_back(LLR.back() < 0);
	}
	else {
		// use freezed value
		decoded_eu.push_back(0);
	}
	LLR.pop_back();

	int ulen = u.size();
	decoded_u.resize(ulen);
	for (int i = 0; i < ulen; ++i) {
		decoded_u[i] = decoded_eu[sorted_channel_ind_ZWi[i].first];
	}
}

void polar::print_decoded_eu() {
	if (decoded_eu.size() != N)	throw "error size of 'decoded_eu' not equals 'N'";
	printf("d_eu : ");
	for (int i = 0; i < N; ++i) {
		printf("%d ", decoded_eu[i]);
	}
	printf("\n");
}

void polar::print_decoded_u() {
	int ulen;
	if (is_CRC16_aided == true)
		ulen = decoded_u.size() - 16;
	else
		ulen = decoded_u.size();
	if (decoded_u.size() != ulen)	throw "error size of 'decoded_u' not equals size of 'u'";
	printf("d_u  : ");
	for (int i = 0; i < ulen; ++i) {
		printf("%d ", decoded_u[i]);
	}
	printf("\n");
}

void polar::SC_decode() {
	decoded_eu.clear();
	LLR.resize(N);
	double scale = 2 / sigma / sigma;
	for (int i = 0; i < N; ++i) {
		LLR[i] = received_v[i] * scale;
	}

	SC_decode_recur(0, N - 1);
	int ulen = u.size();
	decoded_u.resize(ulen);
	for (int i = 0; i < ulen; ++i) {
		decoded_u[i] = decoded_eu[sorted_channel_ind_ZWi[i].first];	// get bit in active channel as decode result
	}
}

void polar::SC_decode_recur(int LLR_start_ind, int LLR_end_ind) {
	int len = LLR_end_ind - LLR_start_ind + 1;
	if (len == 1) {
		if (active_channel_sum[decoded_eu.size() + 1] != active_channel_sum[decoded_eu.size()]) {
			// use LLR value
			decoded_eu.push_back(LLR.back() < 0);
		}
		else {
			// use freezed value
			decoded_eu.push_back(0);
		}
		// compute encoded polar code
		LLR.back() = decoded_eu.back();
		return;
	}

	for (int i = LLR_start_ind; i < LLR_start_ind + len / 2; ++i) {
		// compute LLR with f function for upper part
		LLR.push_back(f(LLR[i], LLR[i + len / 2]));
	}
	// decode upper part of polar code
	SC_decode_recur(LLR_end_ind + 1, LLR.size() - 1);
	for (int i = LLR_start_ind; i < LLR_start_ind + len / 2; ++i) {
		// compute LLR with g function for lower part
		LLR.push_back(g((int)LLR[i + len], LLR[i], LLR[i + len / 2]));
	}
	// decode lower part of polar code
	SC_decode_recur(LLR_end_ind + 1 + len / 2, LLR.size() - 1);
	if (decoded_eu.size() != N) {
		for (int i = LLR_start_ind; i < LLR_start_ind + len / 2; ++i) {
			// compute encoded polar code for decode (upper part)
			LLR[i] = LLR[i + len] != LLR[i + len + len / 2];
			// ... (lower part)
			LLR[i + len / 2] = LLR[i + len / 2 + len];
		}
		for (int i = 0; i < len; ++i) {
			// clear out unused LLR
			LLR.pop_back();
		}
	}
	else {
		// decode done!
		return;
	}
}

void polar::fast_SC_decode() {
	decoded_eu.clear();
	LLR.resize(N);
	decoded_v_tmp.resize(N);
	double scale = 2 / sigma / sigma;
	for (int i = 0; i < N; ++i) {
		LLR[i] = received_v[i] * scale;
	}

	fast_SC_decode_recur(0, N - 1);
	int ulen = u.size();
	decoded_u.resize(ulen);
	for (int i = 0; i < ulen; ++i) {
		decoded_u[i] = decoded_eu[sorted_channel_ind_ZWi[i].first];	// get bit in active channel as decode result
	}
}

void polar::fast_SC_decode_recur(int LLR_start_ind, int LLR_end_ind) {
	int len = LLR_end_ind - LLR_start_ind + 1;
	int start_channel_ind = decoded_eu.size();
	int end_channel_ind = start_channel_ind + len - 1;
	int num_active_channel = active_channel_sum[end_channel_ind + 1] - active_channel_sum[start_channel_ind];
	if (num_active_channel == len) {	// R1 channel
		for (int i = 0; i < len; ++i) {
			// use LLR value
			decoded_v_tmp[i] = (LLR[LLR_start_ind + i] < 0);
			// stroe the encoded polar code for return
			LLR[LLR_start_ind + i] = decoded_v_tmp[i];
		}

		// encode them to decoded_eu
		encode_recur_for_dv(0, len - 1);

		for (int i = 0; i < len; ++i) {
			// store the decode result into decoded_eu
			decoded_eu.push_back(decoded_v_tmp[i]);
		}

		return;
	}
	else if (num_active_channel == len - 1 && active_channel_sum[start_channel_ind + 1] == active_channel_sum[start_channel_ind]) {
		// spc channel
		int xor_decoded_v_tmp = 0;
		for (int i = 0; i < len; ++i) {
			// use LLR value
			decoded_v_tmp[i] = (LLR[LLR_start_ind + i] < 0);

			xor_decoded_v_tmp ^= decoded_v_tmp[i];
		}

		if (xor_decoded_v_tmp == 0) {
			// decoded_v_tmp is ok and don't need to filp
		}
		else {
			// need to filp a bit of decoded_v_tmp to make it valid
			double min_LLR = abs(LLR[LLR_start_ind]);
			int min_LLR_ind = 0;
			for (int i = 1; i < len; ++i) {
				double can = abs(LLR[LLR_start_ind + i]);
				if (min_LLR > can) {
					min_LLR = can;
					min_LLR_ind = i;
				}
			}
			// flip
			decoded_v_tmp[min_LLR_ind] ^= 1;
		}

		for (int i = 0; i < len; ++i) {
			// stroe the encoded polar code for return
			LLR[LLR_start_ind + i] = decoded_v_tmp[i];
		}

		// encode them to decoded_eu
		encode_recur_for_dv(0, len - 1);

		for (int i = 0; i < len; ++i) {
			// store the decode result into decoded_eu
			decoded_eu.push_back(decoded_v_tmp[i]);
		}

		return;
	}/**/
	else if (num_active_channel == 1 && active_channel_sum[end_channel_ind + 1] != active_channel_sum[end_channel_ind]) {
		// rep channel
		double sum_LLR = 0;
		for (int i = 0; i < len; ++i) {
			// compute sum of LLR
			sum_LLR += LLR[LLR_start_ind + i];
		}
		int result_bit = (sum_LLR < 0);
		for (int i = 0; i < len; ++i) {
			LLR[LLR_start_ind + i] = result_bit;
			decoded_eu.push_back(result_bit);	// this decoded_eu is not correct, but decoded_u is correct
		}

		return;
	}
	else if (num_active_channel == 0) {	// R0
		// decoded the bit all into 0
		for (int i = 0; i < len; ++i) {
			LLR[LLR_start_ind + i] = 0;
			decoded_eu.push_back(0);
		}

		return;
	}

	for (int i = LLR_start_ind; i < LLR_start_ind + len / 2; ++i) {
		// compute LLR with f function for upper part
		LLR.push_back(f(LLR[i], LLR[i + len / 2]));
	}
	// decode upper part of polar code
	fast_SC_decode_recur(LLR_end_ind + 1, LLR.size() - 1);
	for (int i = LLR_start_ind; i < LLR_start_ind + len / 2; ++i) {
		// compute LLR with g function for lower part
		LLR.push_back(g((int)LLR[i + len], LLR[i], LLR[i + len / 2]));
	}
	// decode lower part of polar code
	fast_SC_decode_recur(LLR_end_ind + 1 + len / 2, LLR.size() - 1);
	if (decoded_eu.size() != N) {
		for (int i = LLR_start_ind; i < LLR_start_ind + len / 2; ++i) {
			// compute encoded polar code for decode (upper part)
			LLR[i] = LLR[i + len] != LLR[i + len + len / 2];
			// ... (lower part)
			LLR[i + len / 2] = LLR[i + len / 2 + len];
		}
		for (int i = 0; i < len; ++i) {
			// clear out unused LLR
			LLR.pop_back();
		}
	}
	else {
		// decode done!
		return;
	}
}

void polar::encode_recur_for_dv(int start_ind, int end_ind) {
	if (start_ind == end_ind)	return;
	else {
		int len = end_ind - start_ind + 1;
		for (int i = start_ind; i < start_ind + len / 2; ++i) {
			decoded_v_tmp[i] = decoded_v_tmp[i] ^ decoded_v_tmp[i + len / 2];
		}
		encode_recur_for_dv(start_ind, start_ind + len / 2 - 1);
		encode_recur_for_dv(start_ind + len / 2, end_ind);
	}
}

inline double polar::abs(double x) {
	return x > 0 ? x : -x;
}

void polar::set_list_num(int _list_num) {
	list_num = _list_num;
	if (list_num <= 0)	throw "error of 'list_num' less than 0";
	list_num_is_set = true;

	// variables
	//decoded_list_PM;
	//decoded_list_LLR;
	//decoded_list_eu;

}

void polar::SCL_decode() {
	// to be optimized ...
	if (list_num_is_set == false)	throw "error of not seting 'list_num'";

	decoded_list_eu.resize(1);
	decoded_list_eu[0].clear();

	decoded_list_PM.resize(1);
	decoded_list_PM[0] = make_pair(1, 0);

	decoded_list_LLR.resize(1);
	decoded_list_LLR[0].clear();
	//decoded_list_LLR[0].resize(N);
	/*double scale = 2 / sigma / sigma;
	for (int i = 0; i < N; ++i) {
		decoded_list_LLR[0][i] = received_v[i] * scale;
	}/**/
	occupied.resize(list_num);
	is_done.resize(list_num);
	for (int i = 0; i < list_num; ++i) {
		occupied[i] = false;
		is_done[i] = false;
	}

	SCL_decode_recur_father();
	//SCL_decode_recur(0, N - 1);

	if (is_CRC16_aided == false) {
		// find the min PM
		int min_path_ind = decoded_list_PM[0].first - 1;
		double min_PM = decoded_list_PM[0].second;
		for (unsigned int i = 1; i < decoded_list_PM.size(); ++i) {
			if (min_PM > decoded_list_PM[i].second) {
				min_PM = decoded_list_PM[i].second;
				min_path_ind = decoded_list_PM[i].first - 1;
			}
		}
		int ulen = u.size();
		decoded_u.resize(ulen);
		for (int i = 0; i < ulen; ++i) {
			decoded_u[i] = decoded_list_eu[min_path_ind][sorted_channel_ind_ZWi[i].first];	// get bit in active channel as decode result
		}
	}
	else {
		is_SCL_decode_failed = false;
		sort(decoded_list_PM.begin(), decoded_list_PM.end(), \
			[](pair<int, double> p1, pair<int, double> p2) {
				return p1.second < p2.second;
			});
		for (unsigned int j = 0; j < decoded_list_PM.size(); ++j) {
			int ulen = u.size();
			decoded_u.resize(ulen);
			int path_ind = decoded_list_PM[j].first - 1;
			for (int i = 0; i < ulen; ++i) {
				decoded_u[i] = decoded_list_eu[path_ind][sorted_channel_ind_ZWi[i].first];	// get bit in active channel as decode result
			}
			if (is_du_satisfied_crc16()) {
				return;
			}
		}

		// no matching CRC16 in all the list, use the path with minimum path metric ?? or just erase this result
		int ulen = u.size();
		decoded_u.resize(ulen);
		int path_ind = decoded_list_PM[0].first - 1;
		for (int i = 0; i < ulen; ++i) {
			decoded_u[i] = decoded_list_eu[path_ind][sorted_channel_ind_ZWi[i].first];	// get bit in active channel as decode result
		}
		//is_SCL_decode_failed = true;		// erase the result makes error rate 0
	}
}

void polar::SCL_decode_recur_father() {
	LLR.resize(N);
	double scale = 2 / sigma / sigma;
	for (int i = 0; i < N; ++i) {
		LLR[i] = received_v[i] * scale;
	}
	for (int i = 0; i < N / 2; ++i) {
		// compute LLR with f function for upper part
		decoded_list_LLR[0].push_back(f(LLR[i], LLR[i + N / 2]));
		//LLR.push_back(f(LLR[i], LLR[i + len / 2]));
	}
	// decode upper part of polar code
	SCL_decode_recur(0, N / 2 - 1);
	for (unsigned int m = 0; m < decoded_list_LLR.size(); ++m) {		// travel all LLR in list
		for (int i = 0; i < N / 2; ++i) {
			// compute LLR with g function for lower part
			decoded_list_LLR[m][i] = (g((int)decoded_list_LLR[m][i], LLR[i], LLR[i + N / 2]));
		}
		//LLR.push_back(g((int)LLR[i + len], LLR[i], LLR[i + len / 2]));
	}
	// decode lower part of polar code
	SCL_decode_recur(0, N / 2 - 1);

}

void polar::SCL_decode_recur(int LLR_start_ind, int LLR_end_ind) {
	int len = LLR_end_ind - LLR_start_ind + 1;
	if (len == 1) {
		if (active_channel_sum[decoded_list_eu[0].size() + 1] != \
			active_channel_sum[decoded_list_eu[0].size()]) {		// indicate this is an active channel
			// use LLR value
			int s = decoded_list_PM.size();
			for (int i = 0; i < s; ++i) {
				int LLR_ind = decoded_list_PM[i].first;			// just for coding conviniently
				double LLR_back = decoded_list_LLR[LLR_ind - 1].back();
				decoded_list_PM.push_back(make_pair(-LLR_ind, decoded_list_PM[i].second));		// copy LLR, new item's index equals negtive list index - 1
				decoded_list_PM[i].second += LLR_back < 0 ? -LLR_back : 0;		// assume judge to 0
				decoded_list_PM.back().second += LLR_back < 0 ? 0 : LLR_back;	// assume judge to 1, new item
			}
			int ss = decoded_list_PM.size();		// ss should equals 2*s
			if (ss <= list_num) {					// after all case consumed, list number can hold
				// all path kept

				for (int i = s; i < ss; ++i) {
					decoded_list_PM[i].first = i + 1;		// the original index equals list index + 1
					decoded_list_eu.push_back(decoded_list_eu[i - s]);	// copy the corresponding old eu
					decoded_list_eu.back().push_back(1);	// the last eu member is 1 according to new comming assumsion
					decoded_list_LLR.push_back(decoded_list_LLR[i - s]);// copy the corresponding old LLR
					decoded_list_LLR.back().back() = 1;		// change LLR's last value into code value, that is 1
				}

				for (int i = 0; i < s; ++i) {
					decoded_list_eu[i].push_back(0);		// the original eu add one more code 0 
					decoded_list_LLR[i].back() = 0;			// change LLR's last value into code value, that is 0
				}
			}
			else {
				if (s < list_num) {
					// list not full, make it full first
					for (int i = s; i < list_num; ++i) {
						decoded_list_eu.push_back(vector<int>(0));	// copy the corresponding old eu
						decoded_list_LLR.push_back(vector<double>(0));// copy the corresponding old LLR
					}
				}

				// delete some path (without considering complexity), copying the whole LLR
				sort(decoded_list_PM.begin(), decoded_list_PM.end(), \
					[](const pair<int, double>& p1, const pair<int, double>& p2) {
						return p1.second < p2.second;
					});										// sort the path metric to find the minimun list number path metric

				for (int i = list_num; i < ss; ++i) {
					decoded_list_PM.pop_back();				// delet path metric not belong to minimun list number path metric
				}

				for (int f = 0, e = decoded_list_PM.size() - 1; f < e;) {
					while (f < e && decoded_list_PM[f].first > 0) {
						f++;
					}
					while (f < e && decoded_list_PM[e].first < 0) {
						e--;
					}
					if (f < e) {
						pair<int, double> tmp = decoded_list_PM[f];
						decoded_list_PM[f] = decoded_list_PM[e];
						decoded_list_PM[e] = tmp;
						f++;
						e--;
					}
				}

				//sort(decoded_list_PM.begin(), decoded_list_PM.end(), \
				//	[](const pair<int, double>& p1, const pair<int, double>& p2) {
				//		return p1.first > p2.first;
				//	});


				//vector<vector<int>> tmp_decoded_list_eu = decoded_list_eu;
				//vector<vector<double>> tmp_decoded_list_LLR = decoded_list_LLR;
				//for (int i = 0; i < list_num; ++i) {
				//	occupied[i] = false;				// whether the index of list occupied
				//	is_done[i] = false;
				//}
				bool undone_flag = occupied[0];
				bool done_flag = (!undone_flag);
				//occupied.resize(list_num, false);		
				//vector<bool> is_done(list_num, false);
				int not_occupied_starting_ind = 0;

				for (int i = 0; i < list_num; ++i) {		// update the list to have the minimun path metric
					int LLR_ind = decoded_list_PM[i].first;
					if (LLR_ind > 0) {		// original path metric, it's index equals list index + 1
						occupied[LLR_ind - 1] = done_flag;		// old item kept in its same position
						is_done[i] = done_flag;
						//decoded_list_eu[i] = decoded_list_eu[decoded_list_PM[i].first - 1];		// copy corresponding eu
						//decoded_list_eu[i].push_back(0);	// since it is original path metrix, it's decoded result is 0
						//decoded_list_LLR[i] = decoded_list_LLR[decoded_list_PM[i].first - 1];	// copy corresponding LLR, this is time consuming and unnecessary
						//decoded_list_LLR[i].back() = 0;		// add 0 at eu's back for this bit decode result
						//decoded_list_PM[i].first = i + 1;	// update list index value for next update, that is list index + 1 
					}
					else if (occupied[-LLR_ind - 1] == undone_flag) {
						// new path metric, it's index equals negtive list index - 1 and its original index is not occupied in list

						// replace the original list can save copying

						occupied[-LLR_ind - 1] = done_flag;
						is_done[i] = done_flag;

					}
				}

				//for (int i = 0; i < list_num; ++i) {
				//	int LLR_ind = decoded_list_PM[i].first;

				//}

				for (int i = 0; i < list_num; ++i) {
					int LLR_ind = decoded_list_PM[i].first;
					if (!(LLR_ind < 0 && is_done[i] == undone_flag)) {}
					else {			// new path metric, it's index equals negtive list index - 1
						// find a not occupied index in list to replace for new item
						int j;
						for (j = not_occupied_starting_ind; j < list_num; ++j) {
							if (occupied[j] == undone_flag)	break;
						}
						occupied[j] = done_flag;
						is_done[i] = done_flag;
						not_occupied_starting_ind = j + 1;


						// migrating
						decoded_list_eu[j] = decoded_list_eu[-LLR_ind - 1];	// copy corresponding eu
						//decoded_list_eu[i].push_back(1);	// since it is new path metrix, it's decoded result is 1
						decoded_list_LLR[j] = decoded_list_LLR[-LLR_ind - 1];	// copy corresponding LLR
						//decoded_list_LLR[i].back() = 1;		// add 1 at eu's back for this bit decode result
						decoded_list_PM[i].first = -j - 1;		// update list index value for next update, that is list index + 1 

					}
				}


				// fill decoded_list_eu[i].push_back(), decoded_list_LLR[i].back(), decoded_list_PM[i].first
				for (int i = 0; i < list_num; ++i) {
					int LLR_ind = decoded_list_PM[i].first;
					if (LLR_ind > 0) {
						decoded_list_eu[LLR_ind - 1].push_back(0);
						decoded_list_LLR[LLR_ind - 1].back() = 0;
						//decoded_list_PM[i].first = i + 1;
					}
					else {
						decoded_list_PM[i].first = -LLR_ind;
						decoded_list_eu[-LLR_ind - 1].push_back(1);
						decoded_list_LLR[-LLR_ind - 1].back() = 1;
					}
				}

			}
			//decoded_eu.push_back(LLR.back() < 0);
		}
		else {										// this is a freezed channel
			// use freezed value
			int s = decoded_list_PM.size();
			for (int i = 0; i < s; ++i) {
				int LLR_ind = decoded_list_PM[i].first;
				double LLR_back = decoded_list_LLR[LLR_ind - 1].back();
				decoded_list_PM[i].second += LLR_back < 0 ? -LLR_back : 0;
				decoded_list_eu[LLR_ind - 1].push_back(0);		// decode this bit as 0 and add it to the back of eu
				decoded_list_LLR[LLR_ind - 1].back() = 0;		// LLR's back change to decoded bit 0
			}
			//decoded_eu.push_back(0);
		}
		// compute encoded polar code
		//LLR.back() = decoded_eu.back();

		return;
	}
	for (unsigned int m = 0; m < decoded_list_LLR.size(); ++m) {		// travel all LLR in list
		for (int i = LLR_start_ind; i < LLR_start_ind + len / 2; ++i) {
			// compute LLR with f function for upper part
			decoded_list_LLR[m].push_back(f(decoded_list_LLR[m][i], decoded_list_LLR[m][i + len / 2]));
		}
		//LLR.push_back(f(LLR[i], LLR[i + len / 2]));
	}
	// decode upper part of polar code
	SCL_decode_recur(LLR_end_ind + 1, decoded_list_LLR[0].size() - 1);

	for (unsigned int m = 0; m < decoded_list_LLR.size(); ++m) {		// travel all LLR in list
		for (int i = LLR_start_ind; i < LLR_start_ind + len / 2; ++i) {
			// compute LLR with g function for lower part
			decoded_list_LLR[m].push_back(g((int)decoded_list_LLR[m][i + len], \
				decoded_list_LLR[m][i], decoded_list_LLR[m][i + len / 2]));
		}
		//LLR.push_back(g((int)LLR[i + len], LLR[i], LLR[i + len / 2]));
	}

	// decode lower part of polar code
	SCL_decode_recur(LLR_end_ind + 1 + len / 2, decoded_list_LLR[0].size() - 1);

	if (decoded_list_eu[0].size() != N) {								// whether not finished decode
		for (unsigned int m = 0; m < decoded_list_LLR.size(); ++m) {	// travel all LLR in list
			for (int i = LLR_start_ind; i < LLR_start_ind + len / 2; ++i) {
				// compute encoded polar code for decode (upper part)
				decoded_list_LLR[m][i] = decoded_list_LLR[m][i + len] != decoded_list_LLR[m][i + len + len / 2];
				// ... (lower part)
				decoded_list_LLR[m][i + len / 2] = decoded_list_LLR[m][i + len / 2 + len];
			}

			// compute encoded polar code for decode (upper part)
			//LLR[i] = LLR[i + len] != LLR[i + len + len / 2];
			// ... (lower part)
			//LLR[i + len / 2] = LLR[i + len / 2 + len];
		}

		for (unsigned int m = 0; m < decoded_list_LLR.size(); ++m) {	// travel all LLR in list
			for (int i = 0; i < len; ++i) {
				// clear out unused LLR
				decoded_list_LLR[m].pop_back();
			}
			//LLR.pop_back();
		}
	}
	else {
		// decode done!
		return;
	}
}

void polar::set_rate(double _rate) {
	if (_rate <= 0 || _rate >= 1) throw "rate should between 0 and 1";
	rate = _rate;
	rate_is_set = true;
}

void polar::generate_random_u() {
	if (rate_is_set == false)	throw "rate should be set frist using 'set_rate'";
	unsigned int ulen = u.size();
	u.resize((int)round(N * rate));
	random_device rd;
	default_random_engine e(rd());
	bernoulli_distribution bu;

	int ulen_new = u.size();
	for (int i = 0; i < ulen_new; ++i) {
		u[i] = bu(e);
	}
	if (is_CRC16_aided == true)	u_append_crc16();
	if (u.size() > (unsigned int) N)	throw "invalid 'u' size which greater than 2^n";
	if (ulen == u.size()) {}
	else {
		ulen = u.size();
		active_channel_sum.resize(N + 1, 0);
		for (unsigned int i = 0; i < ulen; ++i) {
			active_channel_sum[sorted_channel_ind_ZWi[i].first + 1] = 1;
		}
		for (int i = 1; i <= N; ++i) {
			active_channel_sum[i] += active_channel_sum[i - 1];
		}
	}
}

void polar::information_bit_count() {
	if (is_CRC16_aided == true)
		total_bit += u.size() - 16;
	else
		total_bit += u.size();
}

void polar::SC_error_count() {
	int ulen;
	if (is_CRC16_aided == true)
		ulen = decoded_u.size() - 16;
	else
		ulen = decoded_u.size();
	for (int i = 0; i < ulen; ++i) {
		total_SC_error_bit += (decoded_u[i] != u[i]);
	}
}

void polar::fast_SC_error_count() {
	int ulen;
	if (is_CRC16_aided == true)
		ulen = decoded_u.size() - 16;
	else
		ulen = decoded_u.size();
	for (int i = 0; i < ulen; ++i) {
		total_fast_SC_error_bit += (decoded_u[i] != u[i]);
	}
}

void polar::SCL_error_count() {
	int ulen;
	if (is_CRC16_aided == true)
		ulen = decoded_u.size() - 16;
	else
		ulen = decoded_u.size();

	if (is_SCL_decode_failed == true) {
		SCL_decode_failed_time++;
		SCL_decode_failed_bit += ulen;
	}
	else {
		switch (list_num) {
		case 2:
			for (int i = 0; i < ulen; ++i) {
				total_SCL2_error_bit += (decoded_u[i] != u[i]);
			}
			break;
		case 4:
			for (int i = 0; i < ulen; ++i) {
				total_SCL4_error_bit += (decoded_u[i] != u[i]);
			}
			break;
		case 8:
			for (int i = 0; i < ulen; ++i) {
				total_SCL8_error_bit += (decoded_u[i] != u[i]);
			}
			break;
		case 16:
			for (int i = 0; i < ulen; ++i) {
				total_SCL16_error_bit += (decoded_u[i] != u[i]);
			}
			break;
		case 32:
			for (int i = 0; i < ulen; ++i) {
				total_SCL32_error_bit += (decoded_u[i] != u[i]);
			}
			break;
		default:
			for (int i = 0; i < ulen; ++i) {
				total_SCL_error_bit += (decoded_u[i] != u[i]);
			}
		}
	}
}

double polar::SC_error_rate() {
	if (total_bit == 0)	throw "error of no bit decoded and counted, \
						use 'error_count()' first after decoding";
	return (double)total_SC_error_bit / (double)total_bit;
}

double polar::fast_SC_error_rate() {
	if (total_bit == 0)	throw "error of no bit decoded and counted, \
						use 'error_count()' first after decoding";
	return (double)total_fast_SC_error_bit / (double)total_bit;
}

double polar::SCL_error_rate() {
	if (total_bit == 0)	throw "error of no bit decoded and counted, \
						use 'error_count()' first after decoding";
	int real_SCL_error_bit;
	switch (list_num) {
	case 2:real_SCL_error_bit = total_SCL2_error_bit; break;
	case 4:real_SCL_error_bit = total_SCL4_error_bit; break;
	case 8:real_SCL_error_bit = total_SCL8_error_bit; break;
	case 16:real_SCL_error_bit = total_SCL16_error_bit; break;
	case 32:real_SCL_error_bit = total_SCL32_error_bit; break;
	default:real_SCL_error_bit = total_SCL_error_bit;
	}
	return (double)real_SCL_error_bit / ((double)total_bit - (double)SCL_decode_failed_bit);
}

void polar::clear_error_count() {
	total_bit = 0;
	total_SC_error_bit = 0;
	total_fast_SC_error_bit = 0;
	total_SCL_error_bit = 0;
	total_SCL2_error_bit = 0;
	total_SCL4_error_bit = 0;
	total_SCL8_error_bit = 0;
	total_SCL16_error_bit = 0;
	total_SCL32_error_bit = 0;
	SCL_decode_failed_time = 0;
	total_fast_SCL_error_bit = 0;
}

void polar::set_is_CRC16_aided(bool _is_CRC16_aided) {
	is_CRC16_aided = _is_CRC16_aided;
}

unsigned short polar::u_crc16() {
	int len = u.size();
	// pack vdata

	unsigned short crc16 = 0x0000;
	unsigned short crc_h8, crc_l8;

	int i = 0;
	while (i < len) {
		unsigned char data = 0;
		if (i < len - 7) {
			for (int j = 0; j < 8; ++j) {
				data |= u[i + j] << j;
			}
			i += 8;
		}
		else {
			for (int j = 0; i + j < len; ++j) {
				data |= u[i + j] << j;
			}
			i = len;
		}
		crc_h8 = (crc16 >> 8);
		crc_l8 = (crc16 << 8);
		crc16 = crc_l8 ^ crc_table[crc_h8 ^ data];
	}

	return crc16;
}

void polar::u_append_crc16() {
	unsigned short crc16 = u_crc16();
	for (int i = 0; i < 16; ++i) {
		u.push_back(crc16 & 1);
		crc16 >>= 1;
	}
}

void polar::u_remove_crc16() {
	for (int i = 0; i < 16; ++i) {
		u.pop_back();
	}
}

unsigned short polar::du_get_crc16() {
	unsigned short c = 0;
	int ulen = decoded_u.size();
	for (int i = ulen - 1; i >= ulen - 16; --i) {
		c <<= 1;
		c |= decoded_u[i];
	}
	return c;
}

void polar::du_remove_crc16() {
	for (int i = 0; i < 16; ++i) {
		decoded_u.pop_back();
	}
}

unsigned short polar::du_crc16() {
	int len = decoded_u.size();
	// pack vdata

	unsigned short crc16 = 0x0000;
	unsigned short crc_h8, crc_l8;

	int i = 0;
	while (i < len) {
		unsigned char data = 0;
		if (i < len - 7) {
			for (int j = 0; j < 8; ++j) {
				data |= decoded_u[i + j] << j;
			}
			i += 8;
		}
		else {
			for (int j = 0; i + j < len; ++j) {
				data |= decoded_u[i + j] << j;
			}
			i = len;
		}
		crc_h8 = (crc16 >> 8);
		crc_l8 = (crc16 << 8);
		crc16 = crc_l8 ^ crc_table[crc_h8 ^ data];
	}

	return crc16;
}

void polar::du_append_us(unsigned short us) {
	unsigned short crc16 = us;
	for (int i = 0; i < 16; ++i) {
		decoded_u.push_back(crc16 & 1);
		crc16 >>= 1;
	}
}

bool polar::is_du_satisfied_crc16() {
	unsigned short cr = du_get_crc16();
	du_remove_crc16();
	unsigned short cp = du_crc16();
	du_append_us(cr);
	return cr == cp;
}

double polar::get_time_key_routine() {
	return time_key_routine;
}

void polar::fast_SCL_decode() {
	// to be optimized ...
	if (list_num_is_set == false)	throw "error of not seting 'list_num'";

	decoded_list_eu.resize(1);
	decoded_list_eu[0].clear();

	decoded_list_PM.resize(1);
	decoded_list_PM[0] = make_pair(1, 0);

	decoded_list_LLR.resize(1);
	decoded_list_LLR[0].clear();
	//decoded_list_LLR[0].resize(N);
	/*double scale = 2 / sigma / sigma;
	for (int i = 0; i < N; ++i) {
		decoded_list_LLR[0][i] = received_v[i] * scale;
	}/**/
	occupied.resize(list_num);
	is_done.resize(list_num);
	for (int i = 0; i < list_num; ++i) {
		occupied[i] = false;
		is_done[i] = false;
	}

	fast_SCL_decode_recur_father();
	//SCL_decode_recur(0, N - 1);

	if (is_CRC16_aided == false) {
		// find the min PM
		int min_path_ind = decoded_list_PM[0].first - 1;
		double min_PM = decoded_list_PM[0].second;
		for (unsigned int i = 1; i < decoded_list_PM.size(); ++i) {
			if (min_PM > decoded_list_PM[i].second) {
				min_PM = decoded_list_PM[i].second;
				min_path_ind = decoded_list_PM[i].first - 1;
			}
		}
		int ulen = u.size();
		decoded_u.resize(ulen);
		for (int i = 0; i < ulen; ++i) {
			decoded_u[i] = decoded_list_eu[min_path_ind][sorted_channel_ind_ZWi[i].first];	// get bit in active channel as decode result
		}
	}
	else {
		is_SCL_decode_failed = false;
		sort(decoded_list_PM.begin(), decoded_list_PM.end(), \
			[](pair<int, double> p1, pair<int, double> p2) {
				return p1.second < p2.second;
			});
		for (unsigned int j = 0; j < decoded_list_PM.size(); ++j) {
			int ulen = u.size();
			decoded_u.resize(ulen);
			int path_ind = decoded_list_PM[j].first - 1;
			for (int i = 0; i < ulen; ++i) {
				decoded_u[i] = decoded_list_eu[path_ind][sorted_channel_ind_ZWi[i].first];	// get bit in active channel as decode result
			}
			if (is_du_satisfied_crc16()) {
				return;
			}
		}

		// no matching CRC16 in all the list, use the path with minimum path metric ?? or just erase this result
		if (use_adaptive_L == false) {
			int ulen = u.size();
			decoded_u.resize(ulen);
			int path_ind = decoded_list_PM[0].first - 1;
			for (int i = 0; i < ulen; ++i) {
				decoded_u[i] = decoded_list_eu[path_ind][sorted_channel_ind_ZWi[i].first];	// get bit in active channel as decode result
			}
		}
		else {
			if (list_num == 32) {
				int ulen = u.size();
				decoded_u.resize(ulen);
				int path_ind = decoded_list_PM[0].first - 1;
				for (int i = 0; i < ulen; ++i) {
					decoded_u[i] = decoded_list_eu[path_ind][sorted_channel_ind_ZWi[i].first];	// get bit in active channel as decode result
				}
			}
			else {
				use_adaptive_L = false;			// meaning NOT the end of adaption for L, just for flag
			}
		}
		//is_SCL_decode_failed = true;		// erase the result makes error rate 0
	}
}

void polar::fast_SCL_decode_recur_father() {
	LLR.resize(N);
	double scale = 2 / sigma / sigma;
	for (int i = 0; i < N; ++i) {
		LLR[i] = received_v[i] * scale;
	}
	for (int i = 0; i < N / 2; ++i) {
		// compute LLR with f function for upper part
		decoded_list_LLR[0].push_back(f(LLR[i], LLR[i + N / 2]));
		//LLR.push_back(f(LLR[i], LLR[i + len / 2]));
	}
	// decode upper part of polar code
	fast_SCL_decode_recur(0, N / 2 - 1);
	for (unsigned int m = 0; m < decoded_list_LLR.size(); ++m) {		// travel all LLR in list
		for (int i = 0; i < N / 2; ++i) {
			// compute LLR with g function for lower part
			decoded_list_LLR[m][i] = (g((int)decoded_list_LLR[m][i], LLR[i], LLR[i + N / 2]));
		}
		//LLR.push_back(g((int)LLR[i + len], LLR[i], LLR[i + len / 2]));
	}
	// decode lower part of polar code
	fast_SCL_decode_recur(0, N / 2 - 1);
}

void polar::fast_SCL_decode_recur(int LLR_start_ind, int LLR_end_ind) {
	int len = LLR_end_ind - LLR_start_ind + 1;
	int start_channel_ind = decoded_list_eu[0].size();
	int end_channel_ind = start_channel_ind + len - 1;
	int num_active_channel = active_channel_sum[end_channel_ind + 1] - active_channel_sum[start_channel_ind];
	if (num_active_channel == 0) {	// R0
		// decoded the bit all into 0
		//for (int i = 0; i < len; ++i) {
		//	LLR[LLR_start_ind + i] = 0;
		//	decoded_eu.push_back(0);
		//}

		//return;

		int s = decoded_list_PM.size();
		for (int i = 0; i < s; ++i) {
			int LLR_ind = decoded_list_PM[i].first;
			for (int k = 0; k < len; ++k) {
				double LLR_val = decoded_list_LLR[LLR_ind - 1][LLR_start_ind + k];
				decoded_list_PM[i].second += LLR_val < 0 ? -LLR_val : 0;
				// log(1 + exp(-decoded_list_LLR[LLR_ind - 1][LLR_start_ind + k]));
				decoded_list_eu[LLR_ind - 1].push_back(0);		// decode this bit as 0 and add it to the back of eu
				decoded_list_LLR[LLR_ind - 1][LLR_start_ind + k] = 0;		// LLR's back change to decoded bit 0
			}
		}
		//decoded_eu.push_back(0);

		return;
	}

	if (num_active_channel == 1 && active_channel_sum[end_channel_ind + 1] != active_channel_sum[end_channel_ind]) {
		// rep channel
		//double sum_LLR = 0;
		//for (int i = 0; i < len; ++i) {
			// compute sum of LLR
		//	sum_LLR += LLR[LLR_start_ind + i];
		//}
		//int result_bit = (sum_LLR < 0);
		//for (int i = 0; i < len; ++i) {
		//	LLR[LLR_start_ind + i] = result_bit;
		//	decoded_eu.push_back(result_bit);	// this decoded_eu is not correct, but decoded_u is correct
		//}

		// this is an active channel, otherwise is done on R0 switch above
		// use LLR value
		int s = decoded_list_PM.size();
		for (int i = 0; i < s; ++i) {
			int LLR_ind = decoded_list_PM[i].first;			// just for coding conviniently
			decoded_list_PM.push_back(make_pair(-LLR_ind, decoded_list_PM[i].second));		// copy LLR, new item's index equals negtive list index - 1
			for (int k = 0; k < len; ++k) {
				double LLR_val = decoded_list_LLR[LLR_ind - 1][LLR_start_ind + k];
				decoded_list_PM[i].second += LLR_val < 0 ? -LLR_val : 0;		// assume judge to 0
					// log(1 + exp(-decoded_list_LLR[LLR_ind - 1][LLR_start_ind + k]));
				decoded_list_PM.back().second += LLR_val < 0 ? 0 : LLR_val;		// assume judge to 1, new item
					// log(1 + exp(decoded_list_LLR[LLR_ind - 1][LLR_start_ind + k]));
			}
		}
		int ss = decoded_list_PM.size();		// ss should equals 2*s
		if (ss <= list_num) {					// after all case consumed, list number can hold
			// all path kept

			for (int i = s; i < ss; ++i) {
				decoded_list_PM[i].first = i + 1;		// the original index equals list index + 1
				decoded_list_eu.push_back(decoded_list_eu[i - s]);	// copy the corresponding old eu
				decoded_list_LLR.push_back(decoded_list_LLR[i - s]);// copy the corresponding old LLR
				for (int k = 0; k < len; ++k) {
					decoded_list_eu[i].push_back(1);	// the last eu member is 1 according to new comming assumsion, this decoded_eu is not correct, but decoded_u is correct
					decoded_list_LLR[i][LLR_start_ind + k] = 1;		// change LLR's last value into code value, that is 1
				}
			}

			for (int i = 0; i < s; ++i) {
				for (int k = 0; k < len; ++k) {
					decoded_list_eu[i].push_back(0);		// the original eu add code 0, this decoded_eu is not correct, but decoded_u is correct
					decoded_list_LLR[i][LLR_start_ind + k] = 0;			// change LLR's last value into code value, that is 0
				}
			}
		}
		else {
			if (s < list_num) {
				// list not full, make it full first
				for (int i = s; i < list_num; ++i) {
					decoded_list_eu.push_back(vector<int>(0));	// copy the corresponding old eu
					decoded_list_LLR.push_back(vector<double>(0));// copy the corresponding old LLR
				}
			}

			// delete some path (without considering complexity), copying the whole LLR
			sort(decoded_list_PM.begin(), decoded_list_PM.end(), \
				[](const pair<int, double>& p1, const pair<int, double>& p2) {
					return p1.second < p2.second;
				});										// sort the path metric to find the minimun list number path metric

			for (int i = list_num; i < ss; ++i) {
				decoded_list_PM.pop_back();				// delet path metric not belong to minimun list number path metric
			}

			for (int f = 0, e = decoded_list_PM.size() - 1; f < e;) {
				while (f < e && decoded_list_PM[f].first > 0) {
					f++;
				}
				while (f < e && decoded_list_PM[e].first < 0) {
					e--;
				}
				if (f < e) {
					pair<int, double> tmp = decoded_list_PM[f];
					decoded_list_PM[f] = decoded_list_PM[e];
					decoded_list_PM[e] = tmp;
					f++;
					e--;
				}
			}

			//sort(decoded_list_PM.begin(), decoded_list_PM.end(), \
			//	[](const pair<int, double>& p1, const pair<int, double>& p2) {
			//		return p1.first > p2.first;
			//	});


			//vector<vector<int>> tmp_decoded_list_eu = decoded_list_eu;
			//vector<vector<double>> tmp_decoded_list_LLR = decoded_list_LLR;
			//for (int i = 0; i < list_num; ++i) {
			//	occupied[i] = false;				// whether the index of list occupied
			//	is_done[i] = false;
			//}
			bool undone_flag = occupied[0];
			bool done_flag = (!undone_flag);
			//occupied.resize(list_num, false);		
			//vector<bool> is_done(list_num, false);
			int not_occupied_starting_ind = 0;

			for (int i = 0; i < list_num; ++i) {		// update the list to have the minimun path metric
				int LLR_ind = decoded_list_PM[i].first;
				if (LLR_ind > 0) {		// original path metric, it's index equals list index + 1
					occupied[LLR_ind - 1] = done_flag;		// old item kept in its same position
					is_done[i] = done_flag;
					//decoded_list_eu[i] = decoded_list_eu[decoded_list_PM[i].first - 1];		// copy corresponding eu
					//decoded_list_eu[i].push_back(0);	// since it is original path metrix, it's decoded result is 0
					//decoded_list_LLR[i] = decoded_list_LLR[decoded_list_PM[i].first - 1];	// copy corresponding LLR, this is time consuming and unnecessary
					//decoded_list_LLR[i].back() = 0;		// add 0 at eu's back for this bit decode result
					//decoded_list_PM[i].first = i + 1;	// update list index value for next update, that is list index + 1 
				}
				else if (occupied[-LLR_ind - 1] == undone_flag) {
					// new path metric, it's index equals negtive list index - 1 and its original index is not occupied in list

					// replace the original list can save copying

					occupied[-LLR_ind - 1] = done_flag;
					is_done[i] = done_flag;

				}
			}

			//for (int i = 0; i < list_num; ++i) {
			//	int LLR_ind = decoded_list_PM[i].first;

			//}

			for (int i = 0; i < list_num; ++i) {
				int LLR_ind = decoded_list_PM[i].first;
				if (!(LLR_ind < 0 && is_done[i] == undone_flag)) {}
				else {			// new path metric, it's index equals negtive list index - 1
					// find a not occupied index in list to replace for new item
					int j;
					for (j = not_occupied_starting_ind; j < list_num; ++j) {
						if (occupied[j] == undone_flag)	break;
					}
					occupied[j] = done_flag;
					is_done[i] = done_flag;
					not_occupied_starting_ind = j + 1;


					// migrating
					decoded_list_eu[j] = decoded_list_eu[-LLR_ind - 1];	// copy corresponding eu
					//decoded_list_eu[i].push_back(1);	// since it is new path metrix, it's decoded result is 1
					decoded_list_LLR[j] = decoded_list_LLR[-LLR_ind - 1];	// copy corresponding LLR
					//decoded_list_LLR[i].back() = 1;		// add 1 at eu's back for this bit decode result
					decoded_list_PM[i].first = -j - 1;		// update list index value for next update, that is list index + 1 

				}
			}


			// fill decoded_list_eu[i].push_back(), decoded_list_LLR[i].back(), decoded_list_PM[i].first
			for (int i = 0; i < list_num; ++i) {
				int LLR_ind = decoded_list_PM[i].first;
				if (LLR_ind > 0) {
					for (int k = 0; k < len; ++k) {
						decoded_list_eu[LLR_ind - 1].push_back(0);	// this decoded_eu is not correct, but decoded_u is correct
						decoded_list_LLR[LLR_ind - 1][LLR_start_ind + k] = 0;
					}
					//decoded_list_PM[i].first = i + 1;
				}
				else {
					decoded_list_PM[i].first = -LLR_ind;
					for (int k = 0; k < len; ++k) {
						decoded_list_eu[-LLR_ind - 1].push_back(1);	// this decoded_eu is not correct, but decoded_u is correct
						decoded_list_LLR[-LLR_ind - 1][LLR_start_ind + k] = 1;
					}
				}
			}

		}
		//decoded_eu.push_back(LLR.back() < 0);
		// compute encoded polar code
		//LLR.back() = decoded_eu.back();

		return;
	}/**/

	for (unsigned int m = 0; m < decoded_list_LLR.size(); ++m) {		// travel all LLR in list
		for (int i = LLR_start_ind; i < LLR_start_ind + len / 2; ++i) {
			// compute LLR with f function for upper part
			decoded_list_LLR[m].push_back(f(decoded_list_LLR[m][i], decoded_list_LLR[m][i + len / 2]));
		}
		//LLR.push_back(f(LLR[i], LLR[i + len / 2]));
	}
	// decode upper part of polar code
	fast_SCL_decode_recur(LLR_end_ind + 1, decoded_list_LLR[0].size() - 1);

	for (unsigned int m = 0; m < decoded_list_LLR.size(); ++m) {		// travel all LLR in list
		for (int i = LLR_start_ind; i < LLR_start_ind + len / 2; ++i) {
			// compute LLR with g function for lower part
			decoded_list_LLR[m].push_back(g((int)decoded_list_LLR[m][i + len], \
				decoded_list_LLR[m][i], decoded_list_LLR[m][i + len / 2]));
		}
		//LLR.push_back(g((int)LLR[i + len], LLR[i], LLR[i + len / 2]));
	}

	// decode lower part of polar code
	fast_SCL_decode_recur(LLR_end_ind + 1 + len / 2, decoded_list_LLR[0].size() - 1);

	if (decoded_list_eu[0].size() != N) {								// whether not finished decode
		for (unsigned int m = 0; m < decoded_list_LLR.size(); ++m) {	// travel all LLR in list
			for (int i = LLR_start_ind; i < LLR_start_ind + len / 2; ++i) {
				// compute encoded polar code for decode (upper part)
				decoded_list_LLR[m][i] = decoded_list_LLR[m][i + len] != decoded_list_LLR[m][i + len + len / 2];
				// ... (lower part)
				decoded_list_LLR[m][i + len / 2] = decoded_list_LLR[m][i + len / 2 + len];
			}

			// compute encoded polar code for decode (upper part)
			//LLR[i] = LLR[i + len] != LLR[i + len + len / 2];
			// ... (lower part)
			//LLR[i + len / 2] = LLR[i + len / 2 + len];
		}

		for (unsigned int m = 0; m < decoded_list_LLR.size(); ++m) {	// travel all LLR in list
			for (int i = 0; i < len; ++i) {
				// clear out unused LLR
				decoded_list_LLR[m].pop_back();
			}
			//LLR.pop_back();
		}
	}
	else {
		// decode done!
		return;
	}
}

void polar::fast_SCL_error_count() {
	int ulen;
	if (is_CRC16_aided == true)
		ulen = decoded_u.size() - 16;
	else
		ulen = decoded_u.size();

	if (is_SCL_decode_failed == true) {
		SCL_decode_failed_time++;
		SCL_decode_failed_bit += ulen;
	}
	else {
		if (use_adaptive_L) {
			for (int i = 0; i < ulen; ++i) {
				total_fast_SCL_error_bit += (decoded_u[i] != u[i]);
			}
			return;
		}
		switch (list_num) {
		case 2:
			for (int i = 0; i < ulen; ++i) {
				total_SCL2_error_bit += (decoded_u[i] != u[i]);
			}
			break;
		case 4:
			for (int i = 0; i < ulen; ++i) {
				total_SCL4_error_bit += (decoded_u[i] != u[i]);
			}
			break;
		case 8:
			for (int i = 0; i < ulen; ++i) {
				total_SCL8_error_bit += (decoded_u[i] != u[i]);
			}
			break;
		case 16:
			for (int i = 0; i < ulen; ++i) {
				total_SCL16_error_bit += (decoded_u[i] != u[i]);
			}
			break;
		case 32:
			for (int i = 0; i < ulen; ++i) {
				total_SCL32_error_bit += (decoded_u[i] != u[i]);
			}
			break;
		default:
			for (int i = 0; i < ulen; ++i) {
				total_fast_SCL_error_bit += (decoded_u[i] != u[i]);
			}
		}
	}
}

double polar::fast_SCL_error_rate() {
	if (total_bit == 0)	throw "error of no bit decoded and counted, \
						use 'error_count()' first after decoding";
	int real_SCL_error_bit;
	if (use_adaptive_L) {
		real_SCL_error_bit = total_fast_SCL_error_bit;
	}
	else {
		switch (list_num) {
		case 2:real_SCL_error_bit = total_SCL2_error_bit; break;
		case 4:real_SCL_error_bit = total_SCL4_error_bit; break;
		case 8:real_SCL_error_bit = total_SCL8_error_bit; break;
		case 16:real_SCL_error_bit = total_SCL16_error_bit; break;
		case 32:real_SCL_error_bit = total_SCL32_error_bit; break;
		default:real_SCL_error_bit = total_fast_SCL_error_bit;
		}
	}
	return (double)real_SCL_error_bit / ((double)total_bit - (double)SCL_decode_failed_bit);
}

void polar::adaptive_fast_SCL_decode() {
	if (!is_CRC16_aided) {
		throw "need CRC16, use 'set_is_CRC16_aided()' set CRC16 first";
	}
	else {
		list_num = 1;
		list_num_is_set = true;
		do {
			use_adaptive_L = true;
			fast_SCL_decode();
			list_num <<= 1;
		} while (!use_adaptive_L);
		// 'use_adaptive_L' is true afrer exit
	}
}


class ploar_test{
public:
	int run_example()
	{
		//count_running_time();
		//test_error_rate();
		//test_error_rate2();

		//test_error_rate_adaptive_fast_SCL();
		test_error_rate2({ 2.5 }, 9, 0.5, { 8 }, { 1000 }, true);

		vector<double> SNR = { 1,1.5,2,2.5,3,3.5 };
		const int n = 9;
		const double rate = 0.5;
		vector<int> list_num = { 4,8,16 };
		vector<int> round = { 1000,5000,10000,50000,250000,1000000 };		// to be fixed
		bool is_CRC16_aided = true;
		test_error_rate2(SNR, n, rate, list_num, round, is_CRC16_aided);

		return 0;
	}

	void count_running_time(int start_n, int end_n, double _rate, int _list_num, int _round) {

		printf("\n<--- time consume verse N in one simulation, rate=%f, list_num=%d --->\n", _rate, _list_num);

		printf("N = [ ");
		vector<polar> p;
		for (int i = start_n; i <= end_n; ++i) {
			p.push_back(polar(i));
			p.back().set_rate(_rate);
			p.back().set_SNR(1);				// has no effect to counting running time
			p.back().set_list_num(_list_num);
			//p.back().set_is_CRC16_aided(true);
			if (i != end_n)
				printf("%d, ", 1 << i);
			else
				printf("%d ", 1 << i);
		}
		printf("]\n");

		int plen = p.size();
		clock_t start, end;
		/*
		printf("generate_random_u = [ ");
		for (int j = 0; j < plen; ++j) {
			start = clock();		// starting time
			for (int i = 0; i < _round; ++i) {
				p[j].generate_random_u();
			}
			end = clock();
			if(j!=plen-1)
				printf("%e, ", double(end - start) / CLOCKS_PER_SEC / _round);
			else
				printf("%e ", double(end - start) / CLOCKS_PER_SEC / _round);
		}
		printf("]\n");

		printf("encode = [ ");
		for (int j = 0; j < plen; ++j) {
			start = clock();		// starting time
			for (int i = 0; i < _round; ++i) {
				p[j].encode();
			}
			end = clock();
			if (j != plen - 1)
				printf("%e, ", double(end - start) / CLOCKS_PER_SEC / _round);
			else
				printf("%e ", double(end - start) / CLOCKS_PER_SEC / _round);
		}
		printf("]\n");

		printf("pass_channel = [ ");
		for (int j = 0; j < plen; ++j) {
			start = clock();		// starting time
			for (int i = 0; i < _round; ++i) {
				p[j].pass_channel();
			}
			end = clock();
			if (j != plen - 1)
				printf("%e, ", double(end - start) / CLOCKS_PER_SEC / _round);
			else
				printf("%e ", double(end - start) / CLOCKS_PER_SEC / _round);
		}
		printf("]\n");

		printf("SC_decode = [ ");
		for (int j = 0; j < plen; ++j) {
			start = clock();		// starting time
			for (int i = 0; i < _round; ++i) {
				p[j].SC_decode();
			}
			end = clock();
			if (j != plen - 1)
				printf("%e, ", double(end - start) / CLOCKS_PER_SEC / _round);
			else
				printf("%e ", double(end - start) / CLOCKS_PER_SEC / _round);
		}
		printf("]\n");

		printf("fast_SC_decode = [ ");
		for (int j = 0; j < plen; ++j) {
			start = clock();		// starting time
			for (int i = 0; i < _round; ++i) {
				p[j].fast_SC_decode();
			}
			end = clock();
			if (j != plen - 1)
				printf("%e, ", double(end - start) / CLOCKS_PER_SEC / _round);
			else
				printf("%e ", double(end - start) / CLOCKS_PER_SEC / _round);
		}
		printf("]\n");
		*/

		// just test SCL deocde
		for (int j = 0; j < plen; ++j) {
			p[j].generate_random_u();
			p[j].encode();
			p[j].pass_channel();
		}
		printf("fast_SCL_decode = [ ");
		for (int j = 0; j < plen; ++j) {
			start = clock();		// starting time
			for (int i = 0; i < _round; ++i) {
				p[j].fast_SCL_decode();
			}
			end = clock();
			if (j != plen - 1)
				//printf("%e(%e), ", double(end - start) / CLOCKS_PER_SEC / _round, p[j].get_time_key_routine() / _round);
				printf("%e, ", double(end - start) / CLOCKS_PER_SEC / _round);
			else
				//printf("%e(%e) ", double(end - start) / CLOCKS_PER_SEC / _round, p[j].get_time_key_routine() / _round);
				printf("%e ", double(end - start) / CLOCKS_PER_SEC / _round);
		}
		printf("]\n");
	}

	void test_error_rate(double start_SNR, double end_SNR, double SNR_step, \
		int n, double rate, vector<int> list_num, int _round, bool is_CRC16_aided) {

		polar q(n);
		q.set_rate(rate);
		q.set_is_CRC16_aided(is_CRC16_aided);
		int qlen = 1 + (int)round((end_SNR - start_SNR) / SNR_step);
		vector<double> fast_SC_error_rate;
		vector<vector<double>> fast_SCL_error_rate(list_num.size());
		for (int i = 0; i < qlen; ++i) {
			q.set_SNR(start_SNR + i * SNR_step);
			for (int j = 0; j < _round; ++j) {
				q.generate_random_u();
				q.encode();

				q.modulation();
				q.add_noise();
				//q.pass_channel();

				q.information_bit_count();

				//q.SC_decode();
				//q.SC_error_count();

				q.fast_SC_decode();
				q.fast_SC_error_count();

				for (unsigned int k = 0; k < list_num.size(); ++k) {
					q.set_list_num(list_num[k]);

					//q.SCL_decode();
					//q.SCL_error_count();

					q.fast_SCL_decode();
					q.fast_SCL_error_count();
				}/**/
			}
			//printf("SC_error_rate		: %e\n", q.SC_error_rate());

			fast_SC_error_rate.push_back(q.fast_SC_error_rate());
			for (unsigned int k = 0; k < list_num.size(); ++k) {
				q.set_list_num(list_num[k]);
				fast_SCL_error_rate[k].push_back(q.fast_SCL_error_rate());
			}
			q.clear_error_count();
		}

		printf("error rate verse SNR, N=%d, rate=%f\n", (1 << n), rate);

		printf("SNR = [ ");
		for (int i = 0; i < qlen - 1; ++i) {
			printf("%f, ", start_SNR + i * SNR_step);
		}
		printf("%f", start_SNR + ((double)qlen - 1) * SNR_step);
		printf(" ]\n");

		printf("fast_SC_error_rate = [ ");
		for (int i = 0; i < qlen - 1; ++i) {
			printf("%e, ", fast_SC_error_rate[i]);
		}
		printf("%e", fast_SC_error_rate[qlen - 1]);
		printf(" ]\n");

		for (unsigned int k = 0; k < list_num.size(); ++k) {
			printf("fast_SCL%d_error_rate = [ ", list_num[k]);
			for (int i = 0; i < qlen - 1; ++i) {
				printf("%e, ", fast_SCL_error_rate[k][i]);
			}
			printf("%e", fast_SCL_error_rate[k][qlen - 1]);
			printf(" ]\n");
		}
		/**/
	}

	void test_error_rate2(vector<double> SNR, int n, double rate, vector<int> list_num, \
		vector<int> _round, bool is_CRC16_aided) {

		polar q(n);
		q.set_rate(rate);
		q.set_is_CRC16_aided(is_CRC16_aided);
		int qlen = SNR.size();
		vector<double> fast_SC_error_rate;
		vector<vector<double>> fast_SCL_error_rate(list_num.size());
		for (int i = 0; i < qlen; ++i) {
			q.set_SNR(SNR[i]);
			if (_round[i] == 0) {
				fast_SC_error_rate.push_back(0);
				for (unsigned int k = 0; k < list_num.size(); ++k) {
					q.set_list_num(list_num[k]);
					fast_SCL_error_rate[k].push_back(0);
				}
				q.clear_error_count();
				continue;
			}
			for (int j = 0; j < _round[i]; ++j) {
				q.generate_random_u();
				q.encode();

				q.modulation();
				q.add_noise();
				//q.pass_channel();

				q.information_bit_count();

				//q.SC_decode();
				//q.SC_error_count();

				q.fast_SC_decode();
				q.fast_SC_error_count();

				for (unsigned int k = 0; k < list_num.size(); ++k) {
					q.set_list_num(list_num[k]);

					//q.SCL_decode();
					//q.SCL_error_count();

					q.fast_SCL_decode();
					q.fast_SCL_error_count();
				}/**/
			}
			//printf("SC_error_rate		: %e\n", q.SC_error_rate());

			fast_SC_error_rate.push_back(q.fast_SC_error_rate());
			for (unsigned int k = 0; k < list_num.size(); ++k) {
				q.set_list_num(list_num[k]);
				fast_SCL_error_rate[k].push_back(q.fast_SCL_error_rate());
			}
			q.clear_error_count();
		}

		printf("error rate verse SNR, N=%d, rate=%f\n", (1 << n), rate);

		printf("SNR = [ ");
		for (int i = 0; i < qlen - 1; ++i) {
			printf("%f, ", SNR[i]);
		}
		printf("%f", SNR.back());
		printf(" ]\n");

		printf("fast_SC_error_rate = [ ");
		for (int i = 0; i < qlen - 1; ++i) {
			printf("%e, ", fast_SC_error_rate[i]);
		}
		printf("%e", fast_SC_error_rate[qlen - 1]);
		printf(" ]\n");

		for (unsigned int k = 0; k < list_num.size(); ++k) {
			printf("fast_SCL%d_error_rate = [ ", list_num[k]);
			for (int i = 0; i < qlen - 1; ++i) {
				printf("%e, ", fast_SCL_error_rate[k][i]);
			}
			printf("%e", fast_SCL_error_rate[k][qlen - 1]);
			printf(" ]\n");
		}
		/**/
	}

	void test_error_rate_adaptive_fast_SCL(vector<double> SNR, int n, double rate, vector<int> _round) {
		polar q(n);
		q.set_rate(rate);
		q.set_is_CRC16_aided(true);
		int qlen = SNR.size();
		vector<double> adaptive_fast_SCL_error_rate;
		for (int i = 0; i < qlen; ++i) {
			q.set_SNR(SNR[i]);
			if (_round[i] == 0) {
				adaptive_fast_SCL_error_rate.push_back(0);
				q.clear_error_count();
				continue;
			}
			for (int j = 0; j < _round[i]; ++j) {
				q.generate_random_u();
				q.encode();

				q.modulation();
				q.add_noise();
				//q.pass_channel();

				q.information_bit_count();

				//q.SC_decode();
				//q.SC_error_count();

				//q.fast_SC_decode();
				//q.fast_SC_error_count();

				q.adaptive_fast_SCL_decode();
				q.fast_SCL_error_count();

			}
			//printf("SC_error_rate		: %e\n", q.SC_error_rate());

			adaptive_fast_SCL_error_rate.push_back(q.fast_SCL_error_rate());
			q.clear_error_count();
		}

		printf("error rate verse SNR, N=%d, rate=%f\n", (1 << n), rate);

		printf("SNR = [ ");
		for (int i = 0; i < qlen - 1; ++i) {
			printf("%f, ", SNR[i]);
		}
		printf("%f", SNR.back());
		printf(" ]\n");

		printf("adaptive_fast_SC_error_rate = [ ");
		for (int i = 0; i < qlen - 1; ++i) {
			printf("%e, ", adaptive_fast_SCL_error_rate[i]);
		}
		printf("%e", adaptive_fast_SCL_error_rate[qlen - 1]);
		printf(" ]\n");

		/**/
	}
};
