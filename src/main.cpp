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
    if (input.empty()) {
        throw std::invalid_argument("Input vector is empty");
    }
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
 * @param level Level of wavelet transformation to be performed.
 * 
 * @return A cv::Mat transformed image
 */
cv::Mat wavelet_transform_2D(const cv::Mat& image, const int& level)
{
    int rows = image.rows;
    int cols = image.cols;
    cv::Mat final_trans(rows, cols, CV_32F);

    if (level < 1){
        return image;
    }

    int l = 1;
    while (l <= level) {
        if (l > 1){
            rows = rows/2;
            cols = cols/2;
        }
        
        cv::Mat row_trans(rows, cols, CV_32F);
        cv::Mat row;
        // Apply the transformation for rows
        for (int i = 0; i < rows; ++i){
            if (l == 1){
                row = image.row(i);
            }
            else{
                row = final_trans.row(i).colRange(0,cols);
            }
            vector<float> rowVector(row.begin<float>(), row.end<float>());

            // Performing 1D transformation
            std::vector<float> approx, detail;
            wavelet_transfrom_1D(rowVector, approx, detail);

            // Store approximation in initial col/2 columns 
            // and details in other col/2 columns of the computed row.
            for (int j = 0; j < approx.size(); ++j) {
                row_trans.at<float>(i, j) = approx[j];
                row_trans.at<float>(i, j + cols / 2) = detail[j];
            }
        }

        // Apply the transformation for cols
        for (int j = 0; j < cols; ++j) {
            cv::Mat col = row_trans.col(j);
            vector<float> colVector(col.begin<float>(), col.end<float>());

            // Performing 1D transformation
            std::vector<float> approx, detail;
            wavelet_transfrom_1D(colVector, approx, detail);

            for (int i = 0; i < approx.size(); ++i) {
                final_trans.at<float>(i, j) = approx[i];
                final_trans.at<float>(i + rows / 2, j) = detail[i];
            }
        }
        l = l + 1; 
    }
    return final_trans;
    
}

/**
 * @brief Performs a 1D wavelet transform on an input vector using HAAR.
 *
 * For a given approximations and coefficients vector, this function reconstructs the signal.
 *
 * @param output Reconstructed signal from approximantion and detailes vector.
 * @param approx A reference to the vector to store approximation coefficients.
 * @param detail A reference to the vector to store detail coefficients.
 */
void inverse_wavelet_transfrom_1D(vector<float>& output, 
                          const vector<float>& approx,
                          const vector<float>& detail)
{
    if (approx.empty() || approx.empty()) {
        throw std::invalid_argument("approx or detail vector is empty");
    }

    // size parameter of approx.
    size_t size = approx.size();
    // Adjust the size of the output signal
    output.resize(size*2);
    // Computing the output signal from approx and details
    for (size_t i = 0; i < size; ++i) {
        output[2 * i] = approx[i] + detail[i];
        output[2 * i + 1] = approx[i] - detail[i];
    }
}


/**
 * @brief Performs a multi-level 2D inverse wavelet transform on an input matrix.
 *
 * This function reconstructs the original image from its wavelet coefficients by
 * iteratively applying the inverse wavelet transform to the rows and columns.
 *
 * @param transformed A reference to the matrix containing the transformed wavelet coefficients.
 * @param level The number of decomposition levels to be reversed. 
 *              It should match the number of forward transform levels applied.
 * @return A cv::Mat object representing the reconstructed image.
 */
cv::Mat inverse_wavelet_transform_2D(cv::Mat& transformed, const int& level)
{
    int init_rows = transformed.rows;
    int init_cols = transformed.cols;
    cv::Mat image(init_rows, init_cols, CV_32F);

    int l = level;

    while (l >=1) {
        int rows = init_rows / pow(2, l - 1);
        int cols = init_cols / pow(2, l - 1);

        cv::Mat intermediate(rows, cols, CV_32F);
    
        // Apply inverse transform to columns
        for (int j = 0; j < cols; ++j) {
            std::vector<float> colApprox(rows / 2), colDetail(rows / 2), colOutput(rows);
            for (int i = 0; i < rows / 2; ++i) {
                colApprox[i] = transformed.at<float>(i, j);
                colDetail[i] = transformed.at<float>(i + rows / 2, j);
            }
            // 1D inverse transform
            inverse_wavelet_transfrom_1D(colOutput, colApprox, colDetail);
            for (int i = 0; i < rows; ++i) {
                intermediate.at<float>(i, j) = colOutput[i];
            }
        }

        // Apply inverse transform to rows
        for (int i = 0; i < rows; ++i) {
            std::vector<float> rowApprox(cols / 2), rowDetail(cols / 2), rowOutput(cols);
            for (int j = 0; j < cols / 2; ++j) {
                rowApprox[j] = intermediate.at<float>(i, j);
                rowDetail[j] = intermediate.at<float>(i, j + cols / 2);
            }
            // 1D inverse transform
            inverse_wavelet_transfrom_1D(rowOutput, rowApprox, rowDetail);

            for (int j = 0; j < cols; ++j) {
                transformed.at<float>(i, j) = rowOutput[j];
            }
        }
        l = l - 1;
    }

    return transformed;

} 

int main()
{
    
    // Input parameters

    // image path
    cv::Mat image = cv::imread("../docs/cameraman.png");
    // Level : Level of wavelet transformation
    // * level shouldn't be negative
    // * level shouldn't allow 2^level > image rows or columns
    int level = -2; 

    // checking if the image is valid or not
    if (image.empty()){
        std::cerr << "Error : could not open the image" << std::endl;
    }
    if (image.rows % 2 != 0 || image.cols % 2 != 0) {
        throw std::runtime_error("Image dimensions must be even.");
    }
    image.convertTo(image, CV_32F, 1.0 / 255.0);
    cv::imshow("Cameraman image", image);

    // Wavelet transformation
    cv::Mat WaveTransform = wavelet_transform_2D(image, level);
    cv::imshow("After transformation", WaveTransform);

    // Inverse wavelet transformation
    cv::Mat recon = inverse_wavelet_transform_2D(WaveTransform, level);
    cv::imshow("Reconstructed image", recon);

    cv::waitKey(0);
    return 0;
}