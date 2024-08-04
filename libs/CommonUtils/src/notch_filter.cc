#include "notch_filter.h"

//Constructor
template <typename T>
NotchFilter<T>::NotchFilter(){
	//
}

template <typename T>
NotchFilter<T>::NotchFilter(const array<T, 3> notch_filter_num, const array<T, 3> notch_filter_den):
notch_filter_num_(notch_filter_num),
notch_filter_den_(notch_filter_den){
	//
}

template <typename T>
NotchFilter<T>::~NotchFilter(){
	//
}

template <typename T>
void NotchFilter<T>::UpdateFilterCoeff(const array<T, 3> notch_filter_num, const array<T, 3> notch_filter_den){
	notch_filter_num_ = notch_filter_num;
	notch_filter_den_ = notch_filter_den;
}

template <typename T>
T NotchFilter<T>::Filter(const T u){
	T y;
	y = notch_filter_num_[0]*u + notch_filter_num_[1]*u_1_ + notch_filter_num_[2]*u_2_ - (notch_filter_den_[1]*y_1_ + notch_filter_den_[2]*y_2_);

	u_2_ = u_1_;
	u_1_ = u;

	y_2_ = y_1_;
	y_1_ = y;

	return y;
}

// Explicit template instantiation
template class NotchFilter<float>;
template class NotchFilter<double>;

