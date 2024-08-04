#ifndef NOTCHFILTER_H
#define NOTCHFILTER_H

#include<unistd.h>
#include<atomic>
#include<mutex>
#include<thread>
#include<pthread.h>
#include<iostream>
#include<vector>

using namespace std;

template <typename T>
class NotchFilter{
	public:
		// Constructors
		NotchFilter();
		NotchFilter(const array<T, 3> notch_filter_num, array<T, 3> const notch_filter_den);

		void UpdateFilterCoeff(const array<T, 3> notch_filter_num, const array<T, 3> notch_filter_den);
		T Filter(const T u);

		~NotchFilter();
	private:
		// Filter numerator and denominator. Coefficients are arranged in the order
		//z^0, z^-1, z^-2 .......
		array<T, 3> notch_filter_num_{1, 0, 0};
		array<T, 3> notch_filter_den_{0, 0, 0};

		// Last two filtered outputs
		T y_1_ = 0;
		T y_2_ = 0;

		//Last two inputs
		T u_1_ = 0;
		T u_2_ = 0;
};

#endif