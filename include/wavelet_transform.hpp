#ifndef WAVELET_TRANSFORM_HPP
#define WAVELET_TRANSFORM_HPP

#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
#include <cmath>

using namespace std;

/**
 * @brief Performs a 1D wavelet transform on an input vector using HAAR.
 *
 * For a given input signal, this function computes 
 * the approximation and detail coefficients.
 *
 * @param input An input vector of floating-point values.
 * @param approx A reference to the vector to store approximation coefficients.
 * @param detail A reference to the vector to store detail coefficients.
 * @param lowpass Lowpass filter 
 * @param higlpass High pass filter 
 * @throws std::invalid_argument if input vector is empty.
 */
void wavelet_transfrom_1D(const vector<float>& input, 
                          vector<float>& approx,
                          vector<float>& detail,
                          const vector<float>& lowpass,
                          const vector<float>& highpass);

/**
 * @brief Performs a 2D wavelet transform on an input image using HAAR.
 *
 * For a given input image, this function computes the approximation,
 * row-wise, column-wise, and diagonal details representation.
 *
 * @param image An input image matrix of floating-point values (CV_32F).
 * @param level Level of wavelet transformation to be performed.
 *              Should be greater than 0.
 * @param lowpass Lowpass filter 
 * @param higlpass High pass filter 
 * @return A cv::Mat transformed image.
 */
cv::Mat wavelet_transform_2D(const cv::Mat& image, const int& level);


#endif // WAVELET_TRANSFORM_HPP
