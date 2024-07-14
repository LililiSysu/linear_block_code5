#pragma once
/*****************************************************************//**
 * \file   polynomial_2v.h
 * \brief  polynomial_2v class for any type of coefficients
 *
 * \author lilili
 * \date   September 2022
 *********************************************************************/
#include"../my_lib/Matrix.h"
#include<map>
#include<unordered_map>
using namespace std;

class mod2_special {
public:
	/**
	 * .
	 *
	 * \param n
	 * \param k
	 * \return {n choose k} % 2, use this can prevent overflow, and make the result correct
	 */
	static int n_choose_k(int n, int k) {
		if (n < max_factorial);
		else {
			max_factorial = 2 * n;
			extend();
		}

		return (factorial_factor_2(n) - factorial_factor_2(k) - factorial_factor_2(n - k)) == 0;
	}

	static int max_factorial;
	static Matrix<int> factorial_factor_2;
	static void extend() {
		// compute number of factor 2 of 0! to (max_factorial -1)! and store
		int orig_size = factorial_factor_2.size();
		factorial_factor_2 = factorial_factor_2.combine_right(Matrix<int>(1, max_factorial - orig_size));
		for (int i = orig_size; i < max_factorial; ++i) {
			factorial_factor_2(i) = factorial_factor_2(i - 1) + num_factor_2(i);
		}
		//cout << "factorial_factor_2" << factorial_factor_2;
	}
	/**
	 * .return the "number of factor 2" of num, e.g., 200=2^(4)*25, hence return 4; 64=2^(6), hence return 6
	 */
	static int num_factor_2(int num) {
		if (num <= 0)	return 0;		// only valid for num>0
		int result = 0;
		while (num % 2 == 0) {
			result++;
			num >>= 1;
		}
		return result;
	}
};

int mod2_special::max_factorial = 1001;
// the factorial_factor_2(i) = factor 2 of (i!)
Matrix<int> mod2_special::factorial_factor_2 = Matrix<int>(1, 1001,{
	0,      0,      1,      1,      3,      3,      4,      4,      7,      7,      8,
	8,     10,     10,     11,     11,     15,     15,     16,     16,     18,     18,   
	19,     19,     22,     22,     23,     23,     25,     25,     26,     26,     31,
	31,     32,     32,     34,     34,     35,     35,     38,     38,     39,     39,   
	41,     41,     42,     42,     46,     46,     47,     47,     49,     49,     50,
	50,     53,     53,     54,     54,     56,     56,     57,     57,     63,     63,   
	64,     64,     66,     66,     67,     67,     70,     70,     71,     71,     73,
	73,     74,     74,     78,     78,     79,     79,     81,     81,     82,     82,   
	85,     85,     86,     86,     88,     88,     89,     89,     94,     94,     95,
	95,     97,     97,     98,     98,    101,    101,    102,    102,    104,    104,  
	105,    105,    109,    109,    110,    110,    112,    112,    113,    113,    116,
	116,    117,    117,    119,    119,    120,    120,    127,    127,    128,    128,  
	130,    130,    131,    131,    134,    134,    135,    135,    137,    137,    138,
	138,    142,    142,    143,    143,    145,    145,    146,    146,    149,    149, 
	150,    150,    152,    152,    153,    153,    158,    158,    159,    159,    161,
	161,    162,    162,    165,    165,    166,    166,    168,    168,    169,    169,  
	173,    173,    174,    174,    176,    176,    177,    177,    180,    180,    181,
	181,    183,    183,    184,    184,    190,    190,    191,    191,    193,    193,  
	194,    194,    197,    197,    198,    198,    200,    200,    201,    201,    205,
	205,    206,    206,    208,    208,    209,    209,    212,    212,    213,    213,  
	215,    215,    216,    216,    221,    221,    222,    222,    224,    224,    225,
	225,    228,    228,    229,    229,    231,    231,    232,    232,    236,    236,   
	237,    237,    239,    239,    240,    240,    243,    243,    244,    244,    246,
	246,    247,    247,    255,    255,    256,    256,    258,    258,    259,    259,  
	262,    262,    263,    263,    265,    265,    266,    266,    270,    270,    271,
	271,    273,    273,    274,    274,    277,    277,    278,    278,    280,    280,  
	281,    281,    286,    286,    287,    287,    289,    289,    290,    290,    293,
	293,    294,    294,    296,    296,    297,    297,    301,    301,    302,    302,   
	304,    304,    305,    305,    308,    308,    309,    309,    311,    311,    312,
	312,    318,    318,    319,    319,    321,    321,    322,    322,    325,    325,  
	326,    326,    328,    328,    329,    329,    333,    333,    334,    334,    336,
	336,    337,    337,    340,    340,    341,    341,    343,    343,    344,    344,  
	349,    349,    350,    350,    352,    352,    353,    353,    356,    356,    357,
	357,    359,    359,    360,    360,    364,    364,    365,    365,    367,    367,  
	368,    368,    371,    371,    372,    372,    374,    374,    375,    375,    382,
	382,    383,    383,    385,    385,    386,    386,    389,    389,    390,    390,   
	392,    392,    393,    393,    397,    397,    398,    398,    400,    400,    401,
	401,    404,    404,    405,    405,    407,    407,    408,    408,    413,    413,  
	414,    414,    416,    416,    417,    417,    420,    420,    421,    421,    423,
	423,    424,    424,    428,    428,    429,    429,    431,    431,    432,    432,  
	435,    435,    436,    436,    438,    438,    439,    439,    445,    445,    446,
	446,    448,    448,    449,    449,    452,    452,    453,    453,    455,    455,  
	456,    456,    460,    460,    461,    461,    463,    463,    464,    464,    467,
	467,    468,    468,    470,    470,    471,    471,    476,    476,    477,    477,   
	479,    479,    480,    480,    483,    483,    484,    484,    486,    486,    487,
	487,    491,    491,    492,    492,    494,    494,    495,    495,    498,    498,  
	499,    499,    501,    501,    502,    502,    511,    511,    512,    512,    514,
	514,    515,    515,    518,    518,    519,    519,    521,    521,    522,    522,  
	526,    526,    527,    527,    529,    529,    530,    530,    533,    533,    534, 
	534,    536,    536,    537,    537,    542,    542,    543,    543,    545,    545,  
	546,    546,    549,    549,    550,    550,    552,    552,    553,    553,    557,  
	557,    558,    558,    560,    560,    561,    561,    564,    564,    565,    565,  
	567,    567,    568,    568,    574,    574,    575,    575,    577,    577,    578, 
	578,    581,    581,    582,    582,    584,    584,    585,    585,    589,    589,  
	590,    590,    592,    592,    593,    593,    596,    596,    597,    597,    599, 
	599,    600,    600,    605,    605,    606,    606,    608,    608,    609,    609,  
	612,    612,    613,    613,    615,    615,    616,    616,    620,    620,    621, 
	621,    623,    623,    624,    624,    627,    627,    628,    628,    630,    630,  
	631,    631,    638,    638,    639,    639,    641,    641,    642,    642,    645, 
	645,    646,    646,    648,    648,    649,    649,    653,    653,    654,    654,  
	656,    656,    657,    657,    660,    660,    661,    661,    663,    663,    664,  
	664,    669,    669,    670,    670,    672,    672,    673,    673,    676,    676,  
	677,    677,    679,    679,    680,    680,    684,    684,    685,    685,    687,  
	687,    688,    688,    691,    691,    692,    692,    694,    694,    695,    695,  
	701,    701,    702,    702,    704,    704,    705,    705,    708,    708,    709,  
	709,    711,    711,    712,    712,    716,    716,    717,    717,    719,    719,  
	720,    720,    723,    723,    724,    724,    726,    726,    727,    727,    732, 
	732,    733,    733,    735,    735,    736,    736,    739,    739,    740,    740,  
	742,    742,    743,    743,    747,    747,    748,    748,    750,    750,    751, 
	751,    754,    754,    755,    755,    757,    757,    758,    758,    766,    766,  
	767,    767,    769,    769,    770,    770,    773,    773,    774,    774,    776, 
	776,    777,    777,    781,    781,    782,    782,    784,    784,    785,    785,  
	788,    788,    789,    789,    791,    791,    792,    792,    797,    797,    798,  
	798,    800,    800,    801,    801,    804,    804,    805,    805,    807,    807,  
	808,    808,    812,    812,    813,    813,    815,    815,    816,    816,    819, 
	819,    820,    820,    822,    822,    823,    823,    829,    829,    830,    830,  
	832,    832,    833,    833,    836,    836,    837,    837,    839,    839,    840,  
	840,    844,    844,    845,    845,    847,    847,    848,    848,    851,    851,  
	852,    852,    854,    854,    855,    855,    860,    860,    861,    861,    863, 
	863,    864,    864,    867,    867,    868,    868,    870,    870,    871,    871,  
	875,    875,    876,    876,    878,    878,    879,    879,    882,    882,    883, 
	883,    885,    885,    886,    886,    893,    893,    894,    894,    896,    896,  
	897,    897,    900,    900,    901,    901,    903,    903,    904,    904,    908,  
	908,    909,    909,    911,    911,    912,    912,    915,    915,    916,    916,  
	918,    918,    919,    919,    924,    924,    925,    925,    927,    927,    928, 
	928,    931,    931,    932,    932,    934,    934,    935,    935,    939,    939,  
	940,    940,    942,    942,    943,    943,    946,    946,    947,    947,    949,  
	949,    950,    950,    956,    956,    957,    957,    959,    959,    960,    960,  
	963,    963,    964,    964,    966,    966,    967,    967,    971,    971,    972,  
	972,    974,    974,    975,    975,    978,    978,    979,    979,    981,    981,  
	982,    982,    987,    987,    988,    988,    990,    990,    991,    991,	999
	});


class ipair {
public:
	int i;		// power of x
	int j;		// power of y
	int weighted_degree;
	ipair(int _i = 0, int _j = 0) :i(_i), j(_j), weighted_degree(i* weighted_i + j * weighted_j) {}

	static int weighted_i;
	static int weighted_j;		// this is for computing weighted degree
	static bool is_lex;


	/* define operator {<<,>>} */
	friend ostream& operator << (ostream& out, const ipair& p) {
		out << "(" << p.i << "," << p.j << ")";
		return out;
	}
	friend istream& operator >> (istream& in, ipair& p) {
		in >> p.i >> p.j;
		p.weighted_degree = p.i * weighted_i + p.j * weighted_j;
		return in;
	}

	friend bool operator == (const ipair& p1, const ipair& p2) {
		return (p1.i == p2.i) && (p1.j == p2.j);
	}
	friend bool operator != (const ipair& p1, const ipair& p2) {
		return (p1.i != p2.i) || (p1.j != p2.j);
	}

	friend bool operator < (const ipair& p1, const ipair& p2) {
		if (p1.weighted_degree < p2.weighted_degree) {
			return true;
		}
		else if (p1.weighted_degree == p2.weighted_degree) {
			if (p1.i == p2.i) {
				return false;
			}
			else {
				bool less = p1.i < p2.i;
				return is_lex ? less : (!less);
			}
		}
		else {
			return false;
		}
	}
	friend bool operator > (const ipair& p1, const ipair& p2) {
		return !(p1 < p2);
	}
	friend bool operator <= (const ipair& p1, const ipair& p2) {
		return p1 == p2 || p1 < p2;
	}
	friend bool operator >= (const ipair& p1, const ipair& p2) {
		return p1 == p2 || p1 > p2;
	}

	static Matrix<ipair> under_top(const ipair& top_p, int sum_constrain) {
		Matrix<ipair> result(1, 0, 'v');
		int j = 0;
		while (true) {
			int i = 0;
			while (true) {
				ipair p_test(i, j);
				if (p_test <= top_p && i + j < sum_constrain) {
					result.push_back(p_test);
					i++;
				}
				else {
					break;
				}
			}
			if (i != 0) {
				j++;
			}
			else {
				break;
			}
		}

		result.sort();
		return result;
	}

	/**
	 * .return the Ind(x^K), under (1,v)-revlex order
	 */
	static int A(int K, int v) {
		int r = K % v;
		return (K * K + K * v - r * (v - r)) / 2 / v;
	}
	/**
	 * . return the Ind(y^L), under the (1,v)-revlex order
	 */
	static int B(int L, int v) {
		return L + (L * v * (1 + L)) / 2;
	}
};

int ipair::weighted_i = 1;
int ipair::weighted_j = 1;
bool ipair::is_lex = true;

template<typename T> class polynomial_2v {
private:
	/**
	 * .coefficient of a T polynomial_2v, such that coeff {1,0,0,1} indicates 1 + 0 y + 0 y^2 + 1 y^3
	 * also valid for GF2e element such that coeff {a^3, a^1 ,0 ,1} indicates a^3 + a^1 y + 0 y^2 + y^3
	 * it should be valid also for int and double type, but we never test it.
	 *
	 * we extend this class to bivariate, i.e., Q(x,y)=sum_{i,j>=0}a_{i,j} x^i y^j
	 * in this case the coeff will be 2 dimensional, such that a_{i,j}=coeff(i,j)
	 *
	 */
	Matrix<T> coeff;

public:

	/* define construction functions */

	// this is for conposition onto own written class Matrix
	polynomial_2v(T element = 0) : coeff(T(element)) {}
	/**
	 * .this construction function give a way to initialize with length
	 *
	 * \param len: length of this polynomial_2v
	 * \param type:'0' are usually used to initialize a polynomial_2v
	 */
	polynomial_2v(int row, int col, char type_0_1_i_b_N = '0') :coeff(Matrix<T>(row, col, type_0_1_i_b_N)) {}
	polynomial_2v(const Matrix<T>& _coeff) :coeff(_coeff) {}
	polynomial_2v(const vector<T>& _coeff) :coeff(_coeff) {}

	// Access the individual elements
	T& operator()(const unsigned& pos_ind) {
		return coeff(pos_ind);
	}
	const T& operator()(const unsigned& pos_ind) const {
		return coeff(pos_ind);
	}
	T& operator()(const unsigned& row_ind, const unsigned& col_ind) {
		return coeff(row_ind, col_ind);
	}
	const T& operator()(const unsigned& row_ind, const unsigned& col_ind) const {
		return coeff(row_ind, col_ind);
	}
	
	inline Matrix<T> get_coeff() const {
		return coeff;
	}
	inline int row() const {
		return coeff.row();
	}
	inline int col() const {
		return coeff.col();
	}
	inline int size() const {
		return coeff.size();
	}
	inline const T* ptr() const {
		return coeff.ptr();
	}
	inline polynomial_2v<T> get_part(int row_start_ind, int column_start_ind, int row_end_ind, int column_end_ind) const {
		return coeff.get_part(row_start_ind, column_start_ind, row_end_ind, column_end_ind);
	}
	inline void set_size(int row, int column) {
		coeff.set_size(row, column);
	}

	/**
	 * . append 0 to the left of coeff, change the class
	 *
	 * \param num: the number of 0s to insert left
	 */
	void shift_right(int num) {
		coeff = Matrix<T>(coeff.row(), num, '0').combine_right(coeff);
	}

	/**
	 * . append 0 to the up of coeff, change the class
	 *
	 * \param num: the number of 0s to insert up
	 */
	void shift_down(int num) {
		coeff = Matrix<T>(num, coeff.col(), '0').combine_down(coeff);
	}

	/**
	 * . remove the zero column and rows at the end of the matrix
	 */
	void simplify() {
		int rm = 0;

		int r = row();
		int c = col();
		for (int i = r - 1; i >= 0; --i) {
			int j = 0;
			for (; j < c; ++j) {
				if (coeff(i, j) == 0);
				else{
					break;
				}
			}
			if (j == c) {
				rm++;
			}
			else {
				break;
			}
		}
		r -= rm;
		rm = 0;
		for (int j = c - 1; j >= 0; --j) {
			int i = 0;
			for (; i < r; ++i) {
				if (coeff(i, j) == 0);
				else {
					break;
				}
			}
			if (i == r) {
				rm++;
			}
			else {
				break;
			}
		}
		if (rm == 0) {
			coeff.resize(r, c, true);
			return;
		}
		int new_c = c - rm;

		int ic = c;
		int inew_c = new_c;
		for (int i = 1; i < r; ++i) {
			for (int j = 0; j < new_c; ++j) {
				coeff(inew_c + j) = coeff(ic + j);
			}
			ic+=c;
			inew_c += new_c;
		}

		coeff.resize(r, new_c, true);

		//Matrix<T> result(r, c);
		//for (int i = 0; i < r; ++i) {
		//	for (int j = 0; j < c; ++j) {
		//		result(i, j) = coeff(i, j);
		//	}
		//}
		//coeff = result;

	}

	/**
	 * .
	 *
	 * \return the derivative of polynomial_2v
	 */
	polynomial_2v<T> get_derivative_y() const {
		int r = coeff.row();
		int c = coeff.col();

		if (c > 1) {
			polynomial_2v<T> result(r, c - 1, '0');	// be careful of empty polynomial_2v
			for (int p = 0; p < r; ++p) {
				for (int i = 0; i < c - 1; ++i) {					
					result(p,i) = (i + 1) * coeff(p,i + 1);
				}
			}
			return result;
		}
		else {
			return polynomial_2v<T>(0, 0, '0');
		}
	}

	/* define operator {=, +, -, *, /} */
	//polynomial_2v<T>& operator = (const polynomial_2v<T>& p) {
	//	if (this != &p) {	// judge if set class for itself
	//		coeff = p.coeff;
	//	}
	//	return *this;
	//}
	//polynomial_2v<T>& operator = (const vector<T>& _coeff) {
	//	coeff = Matrix<T>(_coeff);
	//	return *this;
	//}
	//polynomial_2v<T>& operator = (const Matrix<T>& _coeff) {
	//	coeff = _coeff;
	//	return *this;
	//}
	polynomial_2v<T> operator - () const {
		return -coeff;
	}		// add inverse of coeff is itself

	void operator += (const polynomial_2v<T>& p2) {
		int p2r = p2.row();
		int p2c = p2.col();
		if (row() >= p2r && col() >= p2c) {
			for (int i = 0; i < p2r; ++i) {
				for (int j = 0; j < p2c; ++j) {
					coeff(i, j) += p2(i, j);
				}
			}
		}
		else
			(*this) = (*this) - p2;
	}
	void operator -= (const polynomial_2v<T>& p2) {
		int p2r = p2.row();
		int p2c = p2.col();
		if (row() >= p2r && col() >= p2c) {
			for (int i = 0; i < p2r; ++i) {
				for (int j = 0; j < p2c; ++j) {
					coeff(i, j) -= p2(i, j);
				}
			}
		}
		else
			(*this) = (*this) - p2;
	}
	void operator *= (const polynomial_2v<T>& p2) {// this cannot be written as in-place operation
		(*this) = (*this) * p2;
	}
	void operator *= (const T num) {
		int len = size();
		for (int i = 0; i < len; ++i) {
			coeff(i) *= num;
		}
	}
	void operator /= (const polynomial_2v<T>& p2) {
		(*this) = (*this) / p2;
	}
	void operator /= (const T num) {
		int len = size();
		T one_over_num = 1 / num;
		for (int i = 0; i < len; ++i) {
			coeff(i) *= one_over_num;
		}
	}
	void operator %= (const polynomial_2v<T>& p2) {
		(*this) = (*this) % p2;
	}

	friend polynomial_2v<T> operator + (const polynomial_2v<T>& p1, const polynomial_2v<T>& p2) {

		bool flag_max_is_r1 = p1.row() >= p2.row();
		bool flag_max_is_c1 = p1.col() >= p2.col();

		int rmax = flag_max_is_r1 ? p1.row() : p2.row();
		int cmax = flag_max_is_c1 ? p1.col() : p2.col();
		int rmin = flag_max_is_r1 ? p2.row() : p1.row();
		int cmin = flag_max_is_c1 ? p2.col() : p1.col();

		polynomial_2v<T> result(rmax, cmax, '0');
		for (int i = 0; i < rmin; ++i) {
			for (int j = 0; j < cmin; ++j) {
				result(i, j) = p1(i, j) + p2(i, j);
			}
		}
		for (int i = 0; i < rmin; ++i) {
			for (int j = cmin; j < cmax; ++j) {
				result(i, j) = flag_max_is_c1 ? p1(i, j) : p2(i, j);
			}
		}
		for (int i = rmin; i < rmax; ++i) {
			for (int j = 0; j < cmin; ++j) {
				result(i, j) = flag_max_is_r1 ? p1(i, j) : p2(i, j);
			}
		}
		for (int i = rmin; i < rmax; ++i) {
			for (int j = cmin; j < cmax; ++j) {
				result(i, j) = 0;
			}
		}
		result.simplify();
		return result;
	}
	friend polynomial_2v<T> operator - (const polynomial_2v<T>& p1, const polynomial_2v<T>& p2) {

		bool flag_max_is_r1 = p1.row() >= p2.row();
		bool flag_max_is_c1 = p1.col() >= p2.col();

		int rmax = flag_max_is_r1 ? p1.row() : p2.row();
		int cmax = flag_max_is_c1 ? p1.col() : p2.col();
		int rmin = flag_max_is_r1 ? p2.row() : p1.row();
		int cmin = flag_max_is_c1 ? p2.col() : p1.col();

		polynomial_2v<T> result(rmax, cmax, '0');
		for (int i = 0; i < rmin; ++i) {
			for (int j = 0; j < cmin; ++j) {
				result(i, j) = p1(i, j) - p2(i, j);
			}
		}
		for (int i = 0; i < rmin; ++i) {
			for (int j = cmin; j < cmax; ++j) {
				result(i, j) = flag_max_is_c1 ? p1(i, j) : -p2(i, j);
			}
		}
		for (int i = rmin; i < rmax; ++i) {
			for (int j = 0; j < cmin; ++j) {
				result(i, j) = flag_max_is_r1 ? p1(i, j) : -p2(i, j);
			}
		}
		for (int i = rmin; i < rmax; ++i) {
			for (int j = cmin; j < cmax; ++j) {
				result(i, j) = 0;
			}
		}
		result.simplify();
		return result;
	}
	friend polynomial_2v<T> operator * (const polynomial_2v<T>& p1, const polynomial_2v<T>& p2) {
		int r1 = p1.row();
		int c1 = p1.col();
		int r2 = p2.row();
		int c2 = p2.col();

		polynomial_2v<T> result(r1 + r2 - 1, c1 + c2 - 1, '0');
		for (int i = 0; i < r1; ++i) {
			for (int j = 0; j < c1; ++j) {
				for (int t = 0; t < r2; ++t) {
					for (int s = 0; s < c2; ++s) {
						result(i + t, j + s) += p1(i, j) * p2(t, s);
					}
				}
			}
		}
		result.simplify();
		return result;
	}

	friend polynomial_2v<T> operator * (const T num, const polynomial_2v<T>& p) {
		if (num == 0) {
			return polynomial_2v<T>(0, 0, '0');
		}
		int r = p.row();
		int c = p.col();
		int sss = r * c;
		polynomial_2v<T> result(r, c, '0');
		for (int i = 0; i < sss; ++i) {
			result(i) = num * p(i);
		}
		return result;
	}
	friend polynomial_2v<T> operator * (const polynomial_2v<T>& p, const T num) {
		return num * p;
	}
	  
	friend polynomial_2v<T> operator / (polynomial_2v<T> p1, const polynomial_2v<T>& p2) {
		// we assume that p1 and p2 are simplified during calculation
		// if not sure, manully call simplify over p1 and p2, before this function

		// return the quotient only, p1/p2, note that result is not unique, maybe unuseful
		int r1 = p1.row();
		int c1 = p1.col();
		int r2 = p2.row();
		int c2 = p2.col();

		if (r1 >= r2 && c1 >= c2) {

			// find leading monomial of polynomial, by default it will be monomial with max power of x
			int s = c2 - 1;
			for (; s >= 0; s--) {
				if (p2(r2 - 1, s) == 0);
				else {
					break;
				}
			}
			// by now leading polynomial will be in position (r2-1,s)
			polynomial_2v<T> result(r1 - r2 + 1, c1 - c2 + 1, '0');

			for (int i = r1 - 1; i >= r2 - 1; --i) {
				int j = c1 - 1;
				for (; j > c1 - (c2 - s); --j) {
					if (p1(i, j) == 0);
					else {
						// cannot divide any more, return the quotient
						result.simplify();
						return result;
					}
				}

				for (; j >= c2 - 1; --j) {
					if (p1(i, j) == 0);		// skip if the position to eleminate is 0
					else {
						T ratio = p1(i, j) / p2(r2 - 1, s);
						result(i - r2 + 1, j - s) = ratio;
						int p1t = i - r2 + 1;
						int p1s = j - c2 + 1;
						for (int t = 0; t < r2; t++) {
							for (int s = 0; s < c2; s++) {
								p1(p1t + t, p1s + s) -= ratio * p2(t, s);
							}
						}
					}
					
				}
			}
			result.simplify();
			return result;
		}
		else {
			return polynomial_2v<T>(0, 0, '0');
		}
	}
	friend polynomial_2v<T> operator / (const polynomial_2v<T>& p, const T num) {
		int r = p.row();
		int c = p.col();

		int len = r * c;
		polynomial_2v<T> result_coeff(r, c);
		for (int i = 0; i < len; ++i) {
			result_coeff(i) = (1 / num) * p(i);
		}
		return result_coeff;
	}

	friend polynomial_2v<T> operator % (polynomial_2v<T> p1, const polynomial_2v<T>& p2) {
		// we assume that p1 and p2 are simplified during calculation
		// if not sure, manully call simplify over p1 and p2, before this function

		// return residual, p1%p2, note that result is not unique, maybe unuseful
		int r1 = p1.row();
		int c1 = p1.col();
		int r2 = p2.row();
		int c2 = p2.col();

		if (r1 >= r2 && c1 >= c2) {

			// find leading monomial of polynomial, by default it will be monomial with max power of x
			int s = c2 - 1;
			for (; s >= 0; s--) {
				if (p2(r2 - 1, s) == 0);
				else {
					break;
				}
			}
			// by now leading polynomial will be in position (r2-1,s)
			polynomial_2v<T> result(r1 - r2 + 1, c1 - c2 + 1, '0');

			for (int i = r1 - 1; i >= r2 - 1; --i) {
				int j = c1 - 1;
				for (; j > c1 - (c2 - s); --j) {
					if (p1(i, j) == 0);
					else {
						// cannot divide any more, return the quotient
						p1.simplify();
						return p1;
					}
				}

				for (; j >= c2 - 1; --j) {
					if (p1(i, j) == 0);		// skip if the position to eleminate is 0
					else {
						T ratio = p1(i, j) / p2(r2 - 1, s);
						result(i - r2 + 1, j - s) = ratio;
						int p1t = i - r2 + 1;
						int p1s = j - c2 + 1;
						for (int t = 0; t < r2; t++) {
							for (int s = 0; s < c2; s++) {
								p1(p1t + t, p1s + s) -= ratio * p2(t, s);
							}
						}
					}

				}
			}
			p1.simplify();
			return p1;
		}
		else {
			p1.simplify();
			return p1;
		}
	}

	/* define operator {<<,>>} */
	friend ostream& operator << (ostream& out, const polynomial_2v<T>& p) {
		if (p.size() == 0) {		// this should not happen
			T temp1 = { 0 };
			T temp2 = { 1 };
			out << endl << setprecision(4);
			if (temp1 < temp2)
				out << setw(14) << 0 << endl;
			else
				out << setw(6) << 0 << endl;

			return out;
		}
		
		out << p.coeff;
		return out;
	}
	friend istream& operator >> (istream& in, polynomial_2v<T>& p) {
		in >> p.coeff;
		return in;
	}

	/* define operator {==,!=} */
	friend bool operator == (const polynomial_2v<T>& p1, const polynomial_2v<T>& p2) {
		return p1.coeff == p2.coeff;
	}
	friend bool operator != (const polynomial_2v<T>& p1, const polynomial_2v<T>& p2) {
		return p1.coeff != p2.coeff;
	}

	/* define evaluate function */

	/**
	 * . compute the polynomial_2v given x,y in field T
	 *
	 * \param x: the varible x to plug into the polynomial_2v
	 * \param y: the varible y to plug into the polynomial_2v
	 * \return the evaluate result
	 */
	T evaluate(T x, T y) const {
		int r = row();
		int c = col();
		T result = 0;
		T y_pow = 1;
		for (int i = 0; i < r; ++i) {
			T result_i = coeff(i, c - 1);
			for (int j = c - 2; j >= 0; --j) {
				result_i *= y;
				result_i += coeff(i, j);
			}
			result += y_pow * result_i;
			y_pow *= x;
		}	
		return result;
	}

	/* compare is not defined here, return true indicate it is not a ordered structure */
	friend bool operator < (const polynomial_2v<T>& p1, const polynomial_2v<T>& p2) {
		return true;		// should be redefined
	}

	/* 2-dimensional Hasse derivative */
	T Hasse_D(int r, int s, T alpha, T beta) const {
		T result = 0;
		int ro = row();
		int co = col();

		if (r >= ro || s >= co) {
			return result;		// out of polynomial's range, nothing to add
		}

		if (alpha == 0 && beta == 0) {
			return coeff(r, s);
		}
		else if (alpha == 0) {		// 1-dimensional Hasse dirivative, over columns
			for (int j = s; j < co; ++j) {
				result += my::n_choose_k(j, s) * coeff(r, j) * pow(beta, j - s);
			}
		}
		else if (beta == 0) {		// 1-dimensional Hasse dirivative, over rows
			for (int i = r; i < ro; ++i) {
				result += my::n_choose_k(i, r) * coeff(i, s) * pow(alpha, i - r);
			}
		}
		else {
			for (int i = r; i < ro; ++i) {
				for (int j = s; j < co; ++j) {
					result += my::n_choose_k(i, r) * my::n_choose_k(j, s) * coeff(i, j) * pow(alpha, i - r) * pow(beta, j - s);
				}
			}
		}
		return result;
	}

	/* 2-dimemsional Hasse derivative for GF2e */
	T Hasse_D_for_GF2e(int r, int s, T alpha, T beta) const {
		T result = 0;
		int ro = row();
		int co = col();

		if (r >= ro || s >= co) {
			return result;		// out of polynomial's range, nothing to add
		}

		if (alpha == 0 && beta == 0) {
			return coeff(r, s);
		}
		else if (alpha == 0) {		// 1-dimensional Hasse dirivative, over columns
			for (int j = s; j < co; ++j) {
				result += mod2_special::n_choose_k(j, s) * coeff(r, j) * pow(beta, j - s);

			}
		}
		else if (beta == 0) {		// 1-dimensional Hasse dirivative, over rows
			for (int i = r; i < ro; ++i) {
				result += mod2_special::n_choose_k(i, r) * coeff(i, s) * pow(alpha, i - r);

			}
		}
		else {
			for (int i = r; i < ro; ++i) {
				for (int j = s; j < co; ++j) {
					result += mod2_special::n_choose_k(i, r) * mod2_special::n_choose_k(j, s) * coeff(i, j) * pow(alpha, i - r) * pow(beta, j - s);
					/*T new_cp = mod2_special::n_choose_k(i, r)* mod2_special::n_choose_k(j, s)* coeff(i, j)* pow(alpha, i - r)* pow(beta, j - s);
					T old_cp = my::n_choose_k(i, r)* my::n_choose_k(j, s)* coeff(i, j)* pow(alpha, i - r)* pow(beta, j - s);
					if (new_cp != old_cp) {
						cout << "i=" << i << ", j=" << j << ", r=" << r << ", s=" << s << endl;
						cout << "old| new" << endl;
						cout << my::n_choose_k(i, r) <<", " << mod2_special::n_choose_k(i, r) << endl;
						cout << my::n_choose_k(j, s) << ", " << mod2_special::n_choose_k(j, s) << endl;
						cout << "new_cp=" << new_cp << ", old_cp=" << old_cp << endl;
						cout << "-------" << endl;
					}*/
				}
			}
		}
		return result;
	}

	/* find leading nomomial */
	ipair leading_monomial() {
		// first make sure polynomial is simplified 

		int r = row();
		int c = col(); 		
		int max_j = -1;
		ipair leading_monomial(0, 0);

		for (int i = r - 1; i >= 0; --i) {
			for (int j = c - 1; j > max_j; j--) {
				if (coeff(i, j) == 0);
				else {
					ipair candidate_monomial(i, j);

					// should set the right order before doing that
					if (leading_monomial > candidate_monomial);
					else{
						leading_monomial = candidate_monomial;
						max_j = j;
						if (j == c - 1) {
							return leading_monomial;
						}
					}
					break;
				}
			}
		}
		return leading_monomial;
	}

	T evaluate_when_x_eq_0(T y) {
		int c = col();
		T result = coeff(0, c - 1);
		for (int j = c - 2; j >= 0; --j) {
			result *= y;
			result += coeff(0, j);
		}
		return result;
	}

	/**
	 * .for GF, find root by searching all elements in the {0,1,2,...,field_cardinality - 1}
	 */
	Matrix<T> root_when_x_eq_0(int field_cardinality) {
		// we only need to consider the first row
		Matrix<T> result(1, 0, 'v');

		// searching
		for (int i = 0; i < field_cardinality; ++i) {
			T can_y = T(i);
			if (evaluate_when_x_eq_0(can_y) == 0) {
				result.push_back(can_y);
			}
		}
		return result;
	}

	bool is_0_when_y_eq_0() {
		int r = row();
		for (int i = 0; i < r; ++i) {
			if (coeff(i, 0) == 0);
			else {
				return false;
			}
		}
		return true;
	}

	void times_x_plus_alpha(T alpha) {
		int c = coeff.col();
		int r = coeff.row();
		coeff = coeff.combine_down(Matrix<T>(1, c, '0'));	// one more row

		for (int j = 0; j < c; ++j) {
			for (int i = r - 1; i >= 0; --i) {
				coeff(i + 1, j) += coeff(i, j);
				coeff(i, j) *= alpha;
			}

			//auto iter = coeff[j].end();
			//while (iter != coeff[j].begin()) {
			//	iter--;
			//	int i = iter->first;
			//	T store = iter->second;
			//	assign(i + 1, j, (*this)(i + 1, j) + store);		// times x
			//	assign(i, j, store * alpha);						// times alpha
			//}
		}
	}
};

/* following is another interesting implememtation of 2variable polynomial using map, but it has higher complexity, not using now */
struct monomial_order {

	static unsigned weighted_i;
	static unsigned weighted_j;		// this is for computing weighted degree, make sure it >0 unless you know the consequence
	static bool is_lex;

	bool operator() (const pair<int, int>& m1, const pair<int, int>& m2) const {
		int m1_weigthed_degree = m1.first * weighted_i + m1.second * weighted_j;
		int m2_weighted_degree = m2.first * weighted_i + m2.second * weighted_j;

		if (m1_weigthed_degree < m2_weighted_degree) {
			return true;
		}
		else if (m1_weigthed_degree == m2_weighted_degree) {
			if (m1.first == m2.first) {
				return false;
			}
			else {
				bool less = m1.first < m2.first;
				return is_lex ? less : (!less);
			}
		}
		else {
			return false;
		}
	}

	static bool lt(const pair<int, int>& m1, const pair<int, int>& m2) {		// same as operator()
		monomial_order a;
		return a(m1, m2);
	}
};

unsigned monomial_order::weighted_i = 1;
unsigned monomial_order::weighted_j = 1;
bool monomial_order::is_lex = true;

template<typename T> class polynomial_2v_map{
public:	// this is advanced structure in c++
	map<pair<int, int>, T, monomial_order> pmap;

	int max_x_pow;
	int max_y_pow;			// note that this is just a upper bound

	polynomial_2v_map(int _max_x_pow = 0, int _max_y_pow = 0) : max_x_pow(_max_x_pow), max_y_pow(_max_y_pow) {};
	polynomial_2v_map(const Matrix<T>& coeff) {
		int r = coeff.row();
		int c = coeff.col();
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j) {
				assign(i, j, coeff(i, j));
			}
		}
	}

	/**
	 * . to set value in polynomial, you have to use the assign function
	 */
	void assign(const int& x_pow, const int& y_pow, const T& val) {
		if (val == 0) {
			// release the monomial with 0 value, but this will not change max_x_pow and max_y_pow

			// max_x_pow and max_y_pow is a upper limit, which is not the actural max x or y pow!
			
			// (may be we can fix by another 2 maps, indicating that number of elements wiht pow x and pow y, respectively)
			pmap.erase(make_pair(x_pow, y_pow));
		}
		else {

			// update max_x_pow and max_y_pow
			max_x_pow = x_pow < max_x_pow ? max_x_pow : x_pow;
			max_y_pow = y_pow < max_y_pow ? max_y_pow : y_pow;

			// assign value
			pmap[make_pair(x_pow, y_pow)] = val;

		}
	}

	/**
	 * .operator() to fetch elements
	 */
	const T& operator()(const int& x_pow, const int& y_pow) const {
		auto iter = pmap.find(make_pair(x_pow, y_pow));
		if (iter == pmap.end())
			return 0;
		else
			return iter->second;
	}

	friend ostream& operator << (ostream& out, const polynomial_2v_map& p) {
		auto iter = p.pmap.begin();
		out << endl;
		if (iter == p.pmap.end()) {
			out << 0;
		}
		else {
			out << "(" << iter->second << " x^" << (iter->first).first << " y^" << (iter->first).second << ") ";
			++iter;

			while (iter != p.pmap.end()) {
				out << "+(" << iter->second << " x^" << (iter->first).first << " y^" << (iter->first).second << ") "; 
				++iter;
			}
		}		
		cout << endl;
		return out;
	}

	pair<int, int> leading_monomial() const {
		auto iter = pmap.end();
		if (iter != pmap.begin()) {
			iter--;
			return iter->first;
		}
		else {
			return make_pair(0, 0);		// this is for empty polynomial
		}
	}

	T Hasse_D_for_GF2e(int r, int s, T alpha, T beta) const {
		T result = 0;

		if (r > max_x_pow || s > max_y_pow);	// out of polynomial's range, nothing to add
		else if (alpha == 0 && beta == 0) {
			result = (*this)(r, s);
			//return coeff(r, s);
		}
		else if (alpha == 0) {		// 1-dimensional Hasse dirivative, over columns (fix row)
			for (int j = s; j <= max_y_pow; ++j) {		// not suitable for sparse and large matrix, but here in Kotter algorithm, this is okay
				pair<int, int> key_pow = make_pair(r, j);
				auto iter = pmap.find(key_pow);
				if (iter == pmap.end()) {
					if (monomial_order::lt(key_pow, leading_monomial()));
					else {
						break;
					}
				}
				else {
					result += mod2_special::n_choose_k(j, s) * iter->second * pow(beta, j - s);
				}
			}
		}
		else if (beta == 0) {		// 1-dimensional Hasse dirivative, over rows (fix column)
			for (int i = r; i <= max_x_pow; ++i) {		// not suitable for sparse and large matrix, but here in Kotter algorithm, this is okay
				pair<int, int> key_pow = make_pair(i, s);
				auto iter = pmap.find(key_pow);
				if (iter == pmap.end()) {
					if (monomial_order::lt(key_pow, leading_monomial()));
					else {
						break;
					}
				}
				else {
					result += mod2_special::n_choose_k(i, r) * iter->second * pow(alpha, i - r);
				}
			}
		}
		else {
			auto iter = pmap.lower_bound(make_pair(r, s));
			while (iter != pmap.end()) {

				pair<int, int> key_pow = iter->first;

				if (key_pow.first >= r && key_pow.second >= s) {
					result += mod2_special::n_choose_k(key_pow.first, r) * mod2_special::n_choose_k(key_pow.second, s) \
						* iter->second * pow(alpha, key_pow.first - r) * pow(beta, key_pow.second - s);
				}
				iter++;
			}
		}
		return result;
	}

	void operator *= (const T& num) {
		// we assume that for all x != 0, x * y == 0 only if y == 0
		if (num == 0) {
			pmap.clear();
			return;
		}

		for (auto iter = pmap.begin(); iter != pmap.end();++iter) {
			iter->second *= num;
		}
	}

	void operator /= (const T& num) {
		// we assume that for all x != 0, x * y == 0 only if y == 0
		if (num == 0) {
			pmap.clear();
			return;
		}

		for (auto iter = pmap.begin(); iter != pmap.end(); ++iter) {
			iter->second /= num;
		}
	}

	void operator -= (const polynomial_2v_map<T>& p2) {
		for (auto iter = p2.pmap.begin(); iter != p2.pmap.end(); ++iter) {
			int x_pow = (iter->first).first;
			int y_pow = (iter->first).second;
			assign(x_pow, y_pow, (*this)(x_pow, y_pow) - iter->second);
		}
	}

	void operator += (const polynomial_2v_map<T>& p2) {
		for (auto iter = p2.pmap.begin(); iter != p2.pmap.end(); ++iter) {
			int x_pow = (iter->first).first;
			int y_pow = (iter->first).second;
			assign(x_pow, y_pow, (*this)(x_pow, y_pow) + iter->second);
		}
	}

	void times_x_plus_alpha(T alpha) {
		auto iter = pmap.end();
		while (iter != pmap.begin()) {
			iter--;
			pair<int, int> key_pow = iter->first;
			T store = (*this)(key_pow.first, key_pow.second);
			assign(key_pow.first + 1, key_pow.second, (*this)(key_pow.first + 1, key_pow.second) + store);		// times x
			assign(key_pow.first, key_pow.second, store * alpha);		// times alpha
		}
	}

	friend bool operator < (const polynomial_2v_map<T>& p1, const polynomial_2v_map<T>& p2) {
		pair<int, int> p1_lead = p1.leading_monomial();
		pair<int, int> p2_lead = p2.leading_monomial();
		return monomial_order::lt(p1_lead, p2_lead);
	}
};

template<typename T> class polynomial_2v_umap {
public:	// this is advanced structure in c++
	vector<unordered_map<int, T>> coeff;	//coeff[j][i] means coefficient for y^j x^i, for continuous and small y_pow, such as GS decoding
	int max_y_pow;

	// for printing convinent
	polynomial_2v_umap(const Matrix<T>& _coeff) {
		int r = _coeff.row();
		int c = _coeff.col();
		max_y_pow = c - 1;
		coeff = vector<unordered_map<int, T>>(c);
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j) {
				assign(i, j, _coeff(i, j));
			}
		}
	}
	polynomial_2v_umap(): max_y_pow(-1) {}
	void set_max_y_pow(int _max_y_pow) {
		max_y_pow = _max_y_pow;
		coeff = vector<unordered_map<int, T>>(max_y_pow + 1);
	}
	bool contain_Matrix(const Matrix<T>& _coeff) {
		int r = _coeff.row();
		int c = _coeff.col();
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j) {
				if (_coeff(i, j) != (*this)(i, j)) {
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * . to set value in polynomial, you have to use the assign function
	 */
	void assign(const int& x_pow, const int& y_pow, const T& val) {
		// switch (val), to be optimized, need a bool is_zero = val ==0;
		// make sure y_pow is smaller or equal to max_y_pow

		if (val != 0) {
			// assign value
			coeff[y_pow][x_pow] = val;
			
		}
		else {
			// release the monomial with 0 value, but this will not change max_y_pow

			// max_y_pow is a upper limit, which is not the actural max y pow
			coeff[y_pow].erase(x_pow);
		}
	}

	/**
	 * .operator() to fetch elements, usually we donot use it to fetch elements in loop
	 */
	T operator()(const int& x_pow, const int& y_pow) const {
		// make sure y_pow is smaller than max_y_pow
		auto iter = coeff[y_pow].find(x_pow);
		if (iter != coeff[y_pow].end())
			return iter->second;
		else
			return 0;
	}

	friend ostream& operator << (ostream& out, const polynomial_2v_umap& p) {
		cout << endl;
		for (int j = 0; j <= p.max_y_pow; ++j) {
			for (auto iter = p.coeff[j].begin(); iter != p.coeff[j].end(); ++iter) {
				out << "+(" << iter->second << " x^" << iter->first << " y^" << j << ") ";
			}
			cout << endl;
		}
		cout << endl;
		return out;
	}

	pair<int, int> leading_monomial() const {
		pair<int, int> leadm = make_pair(0, 0);
		for (int j = 0; j <= max_y_pow; ++j) {
			pair<int, int> leadj = make_pair(0, 0);
			for (auto iter = coeff[j].begin(); iter != coeff[j].end(); ++iter) {
				pair<int, int> can = make_pair(iter->first, j);
				leadj = leadj.first > can.first ? leadj : can;
			}
			leadm = monomial_order::lt(leadj, leadm) ? leadm : leadj;
		}
		return leadm;
	}

	T Hasse_D_for_GF2e(int r, int s, T alpha, T beta) const {
		T result = 0;

		if (s > max_y_pow);	// out of polynomial's range, nothing to add
		else if (alpha == 0 && beta == 0) {
			result = (*this)(r, s);
			//return coeff(r, s);
		}
		else if (alpha == 0) {		// 1-dimensional Hasse dirivative, over columns (fix row)
			for (int j = s; j <= max_y_pow; ++j) {	// not suitable for sparse and large matrix, but here in Kotter algorithm, this is okay
				auto iter = coeff[j].find(r);
				if (iter != coeff[j].end()) {
					result += mod2_special::n_choose_k(j, s) * iter->second * pow(beta, j - s);
				}
			}
		}
		else if (beta == 0) {		// 1-dimensional Hasse dirivative, over rows (fix column)
			for (auto iter = coeff[s].begin(); iter != coeff[s].end(); ++iter) {
				int i = iter->first;
				if (i >= r) {
					result += mod2_special::n_choose_k(i, r) * iter->second * pow(alpha, i - r);
				}
			}
		}
		else {
			for (int j = s; j <= max_y_pow; ++j) {
				for (auto iter = coeff[j].begin(); iter != coeff[j].end(); ++iter) {
					int i = iter->first;
					if (i >= r) {
						result += mod2_special::n_choose_k(i, r) * mod2_special::n_choose_k(j, s) \
							* iter->second * pow(alpha, i - r) * pow(beta, j - s);
					}
				}
			}
		}
		return result;
	}

	void operator *= (const T& num) {
		// we assume that for all x != 0, x * y == 0 only if y == 0
		if (num == 0) {
			coeff.clear();
			max_y_pow = 0;
			return;
		}

		for (int j = 0; j <= max_y_pow; ++j) {
			for (auto iter = coeff[j].begin(); iter != coeff[j].end(); ++iter) {
				iter->second *= num;
			}
		}
	}

	void operator /= (const T& num) {
		// we assume that for all x != 0, x * y == 0 only if y == 0
		// make sure that num is not zero
		for (int j = 0; j <= max_y_pow; ++j) {
			for (auto iter = coeff[j].begin(); iter != coeff[j].end(); ++iter) {
				iter->second /= num;
			}
		}
	}

	void operator -= (const polynomial_2v_umap<T>& p2) {
		while (p2.max_y_pow > max_y_pow) {
			coeff.push_back(unordered_map<int, T>());
			max_y_pow++;
		}

		for (int j = 0; j <= p2.max_y_pow; ++j) {			
			for (auto iter = p2.coeff[j].begin(); iter != p2.coeff[j].end(); ++iter) {
				int i = iter->first;
				coeff[j][i] -= iter->second;
				if (coeff[j][i] != 0);
				else {
					coeff[j].erase(i);
				}

				//assign(i, j, (*this)(i, j) - iter->second);
				//cout << "(*this)(" << i << ", " << j << ")=" << (*this)(i, j) << endl;
			}
		}
	}

	void operator += (const polynomial_2v_map<T>& p2) {
		while (p2.max_y_pow > max_y_pow) {
			coeff.push_back(unordered_map<int, T>());
			max_y_pow++;
		}

		for (int j = 0; j <= p2.max_y_pow; ++j) {
			for (auto iter = p2.coeff[j].begin(); iter != p2.coeff[j].end(); ++iter) {
				int i = iter->first;
				coeff[j][i] += iter->second;
				if (coeff[j][i] != 0);
				else {
					coeff[j].erase(i);
				}
			}
		}
	}

	void times_x_plus_alpha(T alpha) {

		// coeff*(x+alpha)
		vector<unordered_map<int, T>> coeff2(max_y_pow + 1);
		for (int j = 0; j <= max_y_pow; ++j) {
			for (auto iter = coeff[j].begin(); iter != coeff[j].end(); ++iter) {
				coeff2[j][iter->first + 1] = iter->second;		// times x
			}
		}

		(*this) *= alpha;
		// (*this) += coeff2;
		for (int j = 0; j <= max_y_pow; ++j) {
			for (auto iter = coeff2[j].begin(); iter != coeff2[j].end(); ++iter) {
				int i = iter->first;
				coeff[j][i] += iter->second;
				if (coeff[j][i] != 0);
				else {
					coeff[j].erase(i);
				}
				//assign(i, j, (*this)(i, j) + iter->second);
			}
		}
	}

	friend bool operator < (const polynomial_2v_umap<T>& p1, const polynomial_2v_umap<T>& p2) {
		// time wasting
		pair<int, int> p1_lead = p1.leading_monomial();
		pair<int, int> p2_lead = p2.leading_monomial();
		return monomial_order::lt(p1_lead, p2_lead);
	}
};
