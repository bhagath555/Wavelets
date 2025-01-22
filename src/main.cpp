#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>

using namespace std;



/**
 * @brief Performs a 1D wavelet transform on an input vector using HAAR.
 *
 * For a given input signal, this function computes 
 * the approximation and detail coefficients
 *
 * @param input An input vector of floating-point values.
 * @param approx A reference to the vector to store approximation coefficients.
 * @param detail A reference to the vector to store detail coefficients.
 */
void wavelet_transfrom_1D(const vector<float>& input, 
                          vector<float>& approx,
                          vector<float>& detail)
{

    // size parameter for approx and detail -> half of input.
    size_t size = input.size() / 2;
    
    // Update the size of approx and detail vectors.
    approx.resize(size);
    detail.resize(size);
    // Computing approx and details for given vector
    for (size_t i = 0; i < size; ++i) {
        // Approximation coefficintes 1D
        approx[i] = (input[2*i] + input[2*i + 1] ) * 0.5; 
        // Detail coeffiecients 1D
        detail[i] = (input[2*i] - input[2*i + 1] ) * 0.5; 
    }

}

/**
 * @brief Performs a 2D wavelet transform on an input image using HAAR.
 *
 * For a given input image, this function computes the approximation,
 * row-wise, column-wise, and diagonal details representation.
 *
 * @param image An input vector of floating-point values.
 * @param LL A reference matrix to store approximation coefficients.
 * @param LH A reference matrix to store column-wise details.
 * @param HL A reference matrix to store row-wise details.
 * @param HH A reference matrix to store diagonal details.
 */
void wavelet_transform_2D(const cv::Mat& image, cv::Mat& LL, cv::Mat& LH, cv::Mat& HL, cv::Mat& HH)
{
    int rows = image.rows;
    int cols = image.cols;

    
}

int main()
{
    std::cout << "Testing the template\n";

    cv::Mat image = cv::imread("../docs/cameraman.png");

    if (image.empty()){
        std::cerr << "Error : could not open the image" << std::endl;
    }

    cv::imshow("Cameraman image", image);
    cv::waitKey(0);
    return 0;
}