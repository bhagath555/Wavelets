#ifndef WAVELET_TRANSFORM_HPP
#define WAVELET_TRANSFORM_HPP

#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>

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
 * @throws std::invalid_argument if input vector is empty.
 */
void wavelet_transfrom_1D(const vector<float>& input, 
                          vector<float>& approx,
                          vector<float>& detail);

/**
 * @brief Performs a 2D wavelet transform on an input image using HAAR.
 *
 * For a given input image, this function computes the approximation,
 * row-wise, column-wise, and diagonal details representation.
 *
 * @param image An input image matrix of floating-point values (CV_32F).
 * @param level Level of wavelet transformation to be performed.
 *              Should be greater than 0.
 * @return A cv::Mat transformed image.
 */
cv::Mat wavelet_transform_2D(const cv::Mat& image, const int& level);

/**
 * @brief Performs a 1D inverse wavelet transform on an input vector using HAAR.
 *
 * For given approximation and detail vectors, this function reconstructs the original signal.
 *
 * @param output A reference to the vector to store the reconstructed signal.
 * @param approx A reference to the vector containing approximation coefficients.
 * @param detail A reference to the vector containing detail coefficients.
 * @throws std::invalid_argument if either approx or detail vector is empty.
 */
void inverse_wavelet_transfrom_1D(vector<float>& output, 
                                  const vector<float>& approx,
                                  const vector<float>& detail);

/**
 * @brief Performs a multi-level 2D inverse wavelet transform on an input matrix.
 *
 * This function reconstructs the original image from its wavelet coefficients by
 * iteratively applying the inverse wavelet transform to the rows and columns.
 *
 * @param transformed A reference to the matrix containing the transformed wavelet coefficients.
 * @param level The number of decomposition levels to be reversed. 
 *              Should match the number of forward transform levels applied.
 * @return A cv::Mat object representing the reconstructed image.
 */
cv::Mat inverse_wavelet_transform_2D(cv::Mat& transformed, const int& level);

#endif // WAVELET_TRANSFORM_HPP
