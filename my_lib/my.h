/*****************************************************************//**
 * \file   my.h
 * \brief  class my, this can save interface, easier to extent
 * 
 * \author lilili
 * \date   October 2022
 *********************************************************************/
#pragma once

#include<limits.h>
#include<random>
#include<ctime>
#include<iostream>

// be carefult that we donot using namespace std here

class Complex;

// rand functions of my definition	
class my {
protected:

	static std::uniform_real_distribution<double> u;
	static std::normal_distribution<double> ga;

	static double* qfunc_table;
public:

	// the qfunc as in matlab, has precision of 0.005 over x, and range from -1 to 5
	static double qfunc(double x) {
		int ind = (int)round(my::abs(x) * 200);

		if (ind >= 1000) {
			// reach to limit
			return x > 0 ? 0 : 1;
		}

		if (x > 0) {
			return qfunc_table[ind];
		}
		else {
			return 1- qfunc_table[ind];
		}
	}

	static std::default_random_engine random_engine;

	// a large rand number, ranging from 0 to (1<<30)-1
	static inline int rand_large() {
		// problem: will the random variable be independent during adjacent call of rand()?
		// but advantage: the random seed is determined, every time the program run, we get the
		// same sequence of random variable, good for decoding algorithm comparing
		return (rand() << 15) + rand();		// range: [0,(1<<30)-1]
	}
	// produce a rand number uniformly from (0,1). special: make it not include 0 and 1
	static inline double rand_u() {
		//return ((rand() + 1) / double(RAND_MAX + 2));
		return ((rand_large() + 1.0) / (1073741824/* 1<<30 */ + 1.0));		// make it random enough! much more diversity
	}
	// produce a rand number of normal gaussian distribution, i.e., N(0,1)
	static inline double rand_ga() {
		// if u1=0, this may produce infinity!
		// there is a implicit problem, that is the two adjacent rand_u()
		// may not independent, hence the generated random variable is not true  Gaussian, 
		// you have to call the rand_standard_ga to get the correct random variable
		return sqrt(-2 * log(rand_u())) * cos(2 * my::pi * my::rand_u());
	}
	// random number 0 or 1
	static inline int rand_01() {
		return (rand() % 2);
	}
	// random interger number in range [start, end]
	static inline int rand_int(int start, int end) {
		return start + rand() % (end - start + 1);
	}
	static inline int rand_int() {
		return rand() % 10;
	}

	// advanced rand funtion in c++ 11, recommand to use them, prevent looping a long sequence
	static inline void set_seed_adv(unsigned int _seed = 1) {	// recommand to use time(0)
		random_engine.seed(_seed);
	}
	static inline int rand_large_adv() {
		return random_engine();
	}
	static inline double rand_u_adv() {
		return u(random_engine);
	}
	static inline double rand_ga_adv() {					// return normal distribution random varuable
		return ga(random_engine);
	}
	static inline int rand_01_adv() {
		return random_engine() % 2;
	}
	static inline int rand_int_adv(int start, int end) {	// uniformly distributed random variable in [start, end]
		return start + random_engine() % (end - start + 1);
	}

	// absolute value of my definition
	template<class T> 
	static inline T abs(const T& num) {
		return num >= T(0) ? num : -num;
	}
	// conjugate value of my definition
	template<class T> 
	static inline T conj(const T& num) {
		return num;
	}

	// this implementation should after class complex, special
	static inline Complex abs(const Complex& c);
	static inline Complex conj(const Complex& c);

	/**
	 * .judge if a number is a prime number
	 *
	 * \param _q: number to judge
	 * \return [true] if it is a prime number, [false] if not
	 */
	static bool is_prime(int _q) {
		if (_q <= 1)
			return false;
		if (_q == 2)
			return true;
		if (_q % 2 == 0)
			return false;
		if (_q == 3)
			return true;
		if (_q % 3 == 0)
			return false;
		int s_q = (int)floor(sqrt(_q));
		bool flag = true;
		for (int n = 5; n <= s_q;) {
			if (_q % n != 0) {
				n += flag ? 2 : 4;
				flag ^= flag;
			}
			else {
				return (_q / n == 1);
			}
		}
		return true;
	}

	/**
	 * .
	 *
	 * \param n
	 * \param k
	 * \return number {n choose k}, which is n!/(k!*(n-k)!)
	 */
	static int n_choose_k(int n, int k) {
		if (k <= n - k) {
			// cancell (n-k)!
			int result = 1;
			for (int i = 0; i < k; ++i) {
				result = result * (n - i) / (i + 1);
			}
			return result;
		}
		else {
			// cancell k!, for (n,k) = (n,n-k)			
			return n_choose_k(n, n - k);
		}
	}

	/**
	 * .dealing with number at most 2^64 = 1.8E19
	 *
	 * \param n
	 * \param k
	 * \return number {n choose k}, which is n!/(k!*(n-k)!)
	 */
	static unsigned long long n_choose_k_long(int n, int k) {
		if (k <= n - k) {
			// cancell (n-k)!
			unsigned long long result = 1;
			for (int i = 0; i < k; ++i) {
				result *= unsigned long long(n - i);
				result /= unsigned long long(i + 1);
			}
			return result;
		}
		else {
			// cancell k!, for (n,k) = (n,n-k)			
			return n_choose_k_long(n, n - k);
		}
	}

	/**
	 * .return the lower bound of log2 of num i.e. input 32 (10000) or 35 (10011) return 5
	 */
	static int _log2(int num) {
		int result = -1;
		while (num != 0) {
			num >>= 1;
			result++;
		}
		return result;
	}

	template<class T> 
	static inline T min(const T& x,const T& y) {
		return x < y ? x : y;
	}

	template<class T>
	static inline T max(const T& x,const T& y) {
		return x > y ? x : y;
	}

	static void print_binary(int num, int print_len) {
		int bit_width = _log2(num);
		for (int i = bit_width + 1; i < print_len; ++i) {
			std::cout << 0;
		}
		//std::cout << "bit_width = " << bit_width << std::endl;
		for (int i = bit_width; i >= 0; --i) {
			std::cout << int(((num & (1 << i)) != 0));
		}
	}
	static int rev_bits(int num, int total_bits = 32) {
		int ans = 0;
		total_bits--;
		while (num != 0) {
			ans |= ((num & 1) << total_bits);
			num >>= 1;
			total_bits--;
		}
		return ans;
	}

	static const double pi;
	static const double e;
	static const double zero_approximation;
};

std::default_random_engine my::random_engine;
std::uniform_real_distribution<double> my::u(0,1); 
std::normal_distribution<double> my::ga(0, 1);

const double my::pi = 3.1415926535898;
const double my::e  = 2.7182818284590;
const double my::zero_approximation = 1e-8;

double* my::qfunc_table = new double[1001] {
	0.5,
	4.980E-01,4.960E-01,4.940E-01,4.920E-01,4.900E-01,4.880E-01,4.860E-01,4.840E-01,4.821E-01,4.801E-01,4.781E-01,4.761E-01,4.741E-01,4.721E-01,4.701E-01,4.681E-01,4.661E-01,4.641E-01,4.622E-01,4.602E-01,
	4.582E-01,4.562E-01,4.542E-01,4.522E-01,4.503E-01,4.483E-01,4.463E-01,4.443E-01,4.424E-01,4.404E-01,4.384E-01,4.364E-01,4.345E-01,4.325E-01,4.305E-01,4.286E-01,4.266E-01,4.247E-01,4.227E-01,4.207E-01,
	4.188E-01,4.168E-01,4.149E-01,4.129E-01,4.110E-01,4.090E-01,4.071E-01,4.052E-01,4.032E-01,4.013E-01,3.994E-01,3.974E-01,3.955E-01,3.936E-01,3.917E-01,3.897E-01,3.878E-01,3.859E-01,3.840E-01,3.821E-01,
	3.802E-01,3.783E-01,3.764E-01,3.745E-01,3.726E-01,3.707E-01,3.688E-01,3.669E-01,3.650E-01,3.632E-01,3.613E-01,3.594E-01,3.576E-01,3.557E-01,3.538E-01,3.520E-01,3.501E-01,3.483E-01,3.464E-01,3.446E-01,
	3.427E-01,3.409E-01,3.391E-01,3.372E-01,3.354E-01,3.336E-01,3.318E-01,3.300E-01,3.282E-01,3.264E-01,3.246E-01,3.228E-01,3.210E-01,3.192E-01,3.174E-01,3.156E-01,3.138E-01,3.121E-01,3.103E-01,3.085E-01,
	3.068E-01,3.050E-01,3.033E-01,3.015E-01,2.998E-01,2.981E-01,2.963E-01,2.946E-01,2.929E-01,2.912E-01,2.894E-01,2.877E-01,2.860E-01,2.843E-01,2.826E-01,2.810E-01,2.793E-01,2.776E-01,2.759E-01,2.743E-01,
	2.726E-01,2.709E-01,2.693E-01,2.676E-01,2.660E-01,2.643E-01,2.627E-01,2.611E-01,2.595E-01,2.578E-01,2.562E-01,2.546E-01,2.530E-01,2.514E-01,2.498E-01,2.483E-01,2.467E-01,2.451E-01,2.435E-01,2.420E-01,
	2.404E-01,2.389E-01,2.373E-01,2.358E-01,2.342E-01,2.327E-01,2.312E-01,2.296E-01,2.281E-01,2.266E-01,2.251E-01,2.236E-01,2.221E-01,2.206E-01,2.192E-01,2.177E-01,2.162E-01,2.148E-01,2.133E-01,2.119E-01,
	2.104E-01,2.090E-01,2.075E-01,2.061E-01,2.047E-01,2.033E-01,2.019E-01,2.005E-01,1.991E-01,1.977E-01,1.963E-01,1.949E-01,1.935E-01,1.922E-01,1.908E-01,1.894E-01,1.881E-01,1.867E-01,1.854E-01,1.841E-01,
	1.827E-01,1.814E-01,1.801E-01,1.788E-01,1.775E-01,1.762E-01,1.749E-01,1.736E-01,1.723E-01,1.711E-01,1.698E-01,1.685E-01,1.673E-01,1.660E-01,1.648E-01,1.635E-01,1.623E-01,1.611E-01,1.599E-01,1.587E-01,
	1.574E-01,1.562E-01,1.551E-01,1.539E-01,1.527E-01,1.515E-01,1.503E-01,1.492E-01,1.480E-01,1.469E-01,1.457E-01,1.446E-01,1.434E-01,1.423E-01,1.412E-01,1.401E-01,1.390E-01,1.379E-01,1.368E-01,1.357E-01,
	1.346E-01,1.335E-01,1.324E-01,1.314E-01,1.303E-01,1.292E-01,1.282E-01,1.271E-01,1.261E-01,1.251E-01,1.240E-01,1.230E-01,1.220E-01,1.210E-01,1.200E-01,1.190E-01,1.180E-01,1.170E-01,1.160E-01,1.151E-01,
	1.141E-01,1.131E-01,1.122E-01,1.112E-01,1.103E-01,1.093E-01,1.084E-01,1.075E-01,1.066E-01,1.056E-01,1.047E-01,1.038E-01,1.029E-01,1.020E-01,1.012E-01,1.003E-01,9.940E-02,9.853E-02,9.766E-02,9.680E-02,
	9.595E-02,9.510E-02,9.425E-02,9.342E-02,9.259E-02,9.176E-02,9.094E-02,9.012E-02,8.931E-02,8.851E-02,8.771E-02,8.691E-02,8.613E-02,8.534E-02,8.457E-02,8.379E-02,8.303E-02,8.226E-02,8.151E-02,8.076E-02,
	8.001E-02,7.927E-02,7.853E-02,7.780E-02,7.708E-02,7.636E-02,7.564E-02,7.493E-02,7.423E-02,7.353E-02,7.283E-02,7.215E-02,7.146E-02,7.078E-02,7.011E-02,6.944E-02,6.877E-02,6.811E-02,6.746E-02,6.681E-02,
	6.616E-02,6.552E-02,6.489E-02,6.426E-02,6.363E-02,6.301E-02,6.239E-02,6.178E-02,6.117E-02,6.057E-02,5.997E-02,5.938E-02,5.879E-02,5.821E-02,5.763E-02,5.705E-02,5.648E-02,5.592E-02,5.536E-02,5.480E-02,
	5.425E-02,5.370E-02,5.316E-02,5.262E-02,5.208E-02,5.155E-02,5.102E-02,5.050E-02,4.998E-02,4.947E-02,4.896E-02,4.846E-02,4.796E-02,4.746E-02,4.697E-02,4.648E-02,4.599E-02,4.551E-02,4.504E-02,4.457E-02,
	4.410E-02,4.363E-02,4.317E-02,4.272E-02,4.226E-02,4.182E-02,4.137E-02,4.093E-02,4.049E-02,4.006E-02,3.963E-02,3.920E-02,3.878E-02,3.836E-02,3.795E-02,3.754E-02,3.713E-02,3.673E-02,3.633E-02,3.593E-02,
	3.554E-02,3.515E-02,3.476E-02,3.438E-02,3.400E-02,3.362E-02,3.325E-02,3.288E-02,3.252E-02,3.216E-02,3.180E-02,3.144E-02,3.109E-02,3.074E-02,3.040E-02,3.005E-02,2.971E-02,2.938E-02,2.905E-02,2.872E-02,
	2.839E-02,2.807E-02,2.775E-02,2.743E-02,2.711E-02,2.680E-02,2.650E-02,2.619E-02,2.589E-02,2.559E-02,2.529E-02,2.500E-02,2.471E-02,2.442E-02,2.413E-02,2.385E-02,2.357E-02,2.330E-02,2.302E-02,2.275E-02,
	2.248E-02,2.222E-02,2.195E-02,2.169E-02,2.143E-02,2.118E-02,2.093E-02,2.068E-02,2.043E-02,2.018E-02,1.994E-02,1.970E-02,1.946E-02,1.923E-02,1.899E-02,1.876E-02,1.853E-02,1.831E-02,1.809E-02,1.786E-02,
	1.765E-02,1.743E-02,1.721E-02,1.700E-02,1.679E-02,1.659E-02,1.638E-02,1.618E-02,1.598E-02,1.578E-02,1.558E-02,1.539E-02,1.519E-02,1.500E-02,1.482E-02,1.463E-02,1.444E-02,1.426E-02,1.408E-02,1.390E-02,
	1.373E-02,1.355E-02,1.338E-02,1.321E-02,1.304E-02,1.287E-02,1.271E-02,1.255E-02,1.238E-02,1.222E-02,1.207E-02,1.191E-02,1.176E-02,1.160E-02,1.145E-02,1.130E-02,1.116E-02,1.101E-02,1.087E-02,1.072E-02,
	1.058E-02,1.044E-02,1.031E-02,1.017E-02,1.004E-02,9.903E-03,9.772E-03,9.642E-03,9.514E-03,9.387E-03,9.261E-03,9.137E-03,9.015E-03,8.894E-03,8.774E-03,8.656E-03,8.540E-03,8.424E-03,8.310E-03,8.198E-03,
	8.086E-03,7.976E-03,7.868E-03,7.760E-03,7.654E-03,7.549E-03,7.446E-03,7.344E-03,7.243E-03,7.143E-03,7.044E-03,6.947E-03,6.851E-03,6.756E-03,6.662E-03,6.569E-03,6.478E-03,6.387E-03,6.298E-03,6.210E-03,
	6.123E-03,6.037E-03,5.952E-03,5.868E-03,5.785E-03,5.703E-03,5.622E-03,5.543E-03,5.464E-03,5.386E-03,5.309E-03,5.234E-03,5.159E-03,5.085E-03,5.012E-03,4.940E-03,4.869E-03,4.799E-03,4.730E-03,4.661E-03,
	4.594E-03,4.527E-03,4.461E-03,4.396E-03,4.332E-03,4.269E-03,4.207E-03,4.145E-03,4.085E-03,4.025E-03,3.965E-03,3.907E-03,3.849E-03,3.793E-03,3.736E-03,3.681E-03,3.626E-03,3.573E-03,3.519E-03,3.467E-03,
	3.415E-03,3.364E-03,3.314E-03,3.264E-03,3.215E-03,3.167E-03,3.119E-03,3.072E-03,3.026E-03,2.980E-03,2.935E-03,2.890E-03,2.846E-03,2.803E-03,2.760E-03,2.718E-03,2.676E-03,2.635E-03,2.595E-03,2.555E-03,
	2.516E-03,2.477E-03,2.439E-03,2.401E-03,2.364E-03,2.327E-03,2.291E-03,2.256E-03,2.221E-03,2.186E-03,2.152E-03,2.118E-03,2.085E-03,2.052E-03,2.020E-03,1.988E-03,1.957E-03,1.926E-03,1.896E-03,1.866E-03,
	1.836E-03,1.807E-03,1.778E-03,1.750E-03,1.722E-03,1.695E-03,1.668E-03,1.641E-03,1.615E-03,1.589E-03,1.563E-03,1.538E-03,1.513E-03,1.489E-03,1.465E-03,1.441E-03,1.418E-03,1.395E-03,1.372E-03,1.350E-03,
	1.328E-03,1.306E-03,1.285E-03,1.264E-03,1.243E-03,1.223E-03,1.203E-03,1.183E-03,1.163E-03,1.144E-03,1.125E-03,1.107E-03,1.088E-03,1.070E-03,1.053E-03,1.035E-03,1.018E-03,1.001E-03,9.841E-04,9.676E-04,
	9.514E-04,9.354E-04,9.197E-04,9.043E-04,8.890E-04,8.740E-04,8.593E-04,8.447E-04,8.304E-04,8.164E-04,8.025E-04,7.888E-04,7.754E-04,7.622E-04,7.492E-04,7.364E-04,7.238E-04,7.114E-04,6.992E-04,6.871E-04,
	6.753E-04,6.637E-04,6.522E-04,6.410E-04,6.299E-04,6.190E-04,6.082E-04,5.976E-04,5.873E-04,5.770E-04,5.670E-04,5.571E-04,5.473E-04,5.377E-04,5.283E-04,5.190E-04,5.099E-04,5.009E-04,4.921E-04,4.834E-04,
	4.749E-04,4.665E-04,4.582E-04,4.501E-04,4.421E-04,4.342E-04,4.265E-04,4.189E-04,4.114E-04,4.041E-04,3.968E-04,3.897E-04,3.827E-04,3.758E-04,3.691E-04,3.624E-04,3.559E-04,3.495E-04,3.431E-04,3.369E-04,
	3.308E-04,3.248E-04,3.189E-04,3.131E-04,3.074E-04,3.018E-04,2.963E-04,2.909E-04,2.855E-04,2.803E-04,2.751E-04,2.701E-04,2.651E-04,2.602E-04,2.554E-04,2.507E-04,2.461E-04,2.415E-04,2.370E-04,2.326E-04,
	2.283E-04,2.241E-04,2.199E-04,2.158E-04,2.117E-04,2.078E-04,2.039E-04,2.001E-04,1.963E-04,1.926E-04,1.890E-04,1.854E-04,1.819E-04,1.785E-04,1.751E-04,1.718E-04,1.685E-04,1.653E-04,1.622E-04,1.591E-04,
	1.561E-04,1.531E-04,1.502E-04,1.473E-04,1.445E-04,1.417E-04,1.390E-04,1.363E-04,1.337E-04,1.311E-04,1.286E-04,1.261E-04,1.237E-04,1.213E-04,1.189E-04,1.166E-04,1.144E-04,1.121E-04,1.099E-04,1.078E-04,
	1.057E-04,1.036E-04,1.016E-04,9.961E-05,9.766E-05,9.574E-05,9.386E-05,9.201E-05,9.020E-05,8.842E-05,8.667E-05,8.496E-05,8.327E-05,8.162E-05,8.000E-05,7.841E-05,7.685E-05,7.532E-05,7.382E-05,7.235E-05,
	7.090E-05,6.948E-05,6.809E-05,6.673E-05,6.539E-05,6.407E-05,6.278E-05,6.152E-05,6.028E-05,5.906E-05,5.786E-05,5.669E-05,5.554E-05,5.442E-05,5.331E-05,5.223E-05,5.116E-05,5.012E-05,4.910E-05,4.810E-05,
	4.711E-05,4.615E-05,4.520E-05,4.427E-05,4.336E-05,4.247E-05,4.160E-05,4.074E-05,3.990E-05,3.908E-05,3.827E-05,3.747E-05,3.670E-05,3.594E-05,3.519E-05,3.446E-05,3.374E-05,3.304E-05,3.235E-05,3.167E-05,
	3.101E-05,3.036E-05,2.972E-05,2.910E-05,2.849E-05,2.789E-05,2.730E-05,2.673E-05,2.616E-05,2.561E-05,2.507E-05,2.454E-05,2.402E-05,2.351E-05,2.301E-05,2.252E-05,2.204E-05,2.157E-05,2.111E-05,2.066E-05,
	2.022E-05,1.978E-05,1.936E-05,1.894E-05,1.854E-05,1.814E-05,1.775E-05,1.737E-05,1.699E-05,1.662E-05,1.626E-05,1.591E-05,1.557E-05,1.523E-05,1.490E-05,1.458E-05,1.426E-05,1.395E-05,1.364E-05,1.335E-05,
	1.305E-05,1.277E-05,1.249E-05,1.222E-05,1.195E-05,1.168E-05,1.143E-05,1.118E-05,1.093E-05,1.069E-05,1.045E-05,1.022E-05,9.995E-06,9.774E-06,9.557E-06,9.345E-06,9.137E-06,8.934E-06,8.735E-06,8.540E-06,
	8.349E-06,8.163E-06,7.980E-06,7.801E-06,7.627E-06,7.455E-06,7.288E-06,7.124E-06,6.964E-06,6.807E-06,6.653E-06,6.503E-06,6.356E-06,6.212E-06,6.072E-06,5.934E-06,5.799E-06,5.668E-06,5.539E-06,5.413E-06,
	5.289E-06,5.169E-06,5.050E-06,4.935E-06,4.822E-06,4.712E-06,4.604E-06,4.498E-06,4.395E-06,4.294E-06,4.195E-06,4.098E-06,4.003E-06,3.911E-06,3.821E-06,3.732E-06,3.646E-06,3.561E-06,3.478E-06,3.398E-06,
	3.319E-06,3.241E-06,3.166E-06,3.092E-06,3.020E-06,2.949E-06,2.880E-06,2.813E-06,2.747E-06,2.682E-06,2.619E-06,2.558E-06,2.497E-06,2.439E-06,2.381E-06,2.325E-06,2.270E-06,2.216E-06,2.164E-06,2.112E-06,
	2.062E-06,2.013E-06,1.965E-06,1.919E-06,1.873E-06,1.828E-06,1.785E-06,1.742E-06,1.700E-06,1.660E-06,1.620E-06,1.581E-06,1.543E-06,1.506E-06,1.470E-06,1.434E-06,1.400E-06,1.366E-06,1.333E-06,1.301E-06,
	1.269E-06,1.239E-06,1.209E-06,1.179E-06,1.151E-06,1.123E-06,1.095E-06,1.069E-06,1.043E-06,1.017E-06,9.922E-07,9.680E-07,9.443E-07,9.211E-07,8.985E-07,8.765E-07,8.549E-07,8.339E-07,8.134E-07,7.933E-07,
	7.738E-07,7.547E-07,7.360E-07,7.178E-07,7.000E-07,6.827E-07,6.657E-07,6.492E-07,6.331E-07,6.173E-07,6.019E-07,5.869E-07,5.723E-07,5.580E-07,5.440E-07,5.304E-07,5.171E-07,5.042E-07,4.915E-07,4.792E-07,
	4.671E-07,4.554E-07,4.439E-07,4.327E-07,4.218E-07,4.111E-07,4.008E-07,3.906E-07,3.807E-07,3.711E-07,3.617E-07,3.525E-07,3.435E-07,3.348E-07,3.262E-07,3.179E-07,3.098E-07,3.019E-07,2.942E-07,2.867E-07,
};
