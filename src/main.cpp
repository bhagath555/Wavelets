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
 * , cv::Mat& LL, cv::Mat& LH, cv::Mat& HL, cv::Mat& HH
 */
cv::Mat wavelet_transform_2D(const cv::Mat& image, const int& level)
{
    int rows = image.rows;
    int cols = image.cols;
    cv::Mat final_trans(rows, cols, CV_32F);

    int l = 1;
    while (l <= level) {

        cout << "Level " << l << endl;
        if (l > 1){
            rows = rows/2;
            cols = cols/2;
        }
        
        cv::Mat row_trans(rows, cols, CV_32F);
        cv::Mat row;
        cout << "Row initiation" << endl;
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

    // // Extracting approximation and details
    // LL = final_trans(cv::Rect(0,      0, cols/2, rows /2)).clone();
    // LH = final_trans(cv::Rect(cols/2, 0, cols/2, rows /2)).clone();
    // HL = final_trans(cv::Rect(0, rows/2, cols/2, rows /2)).clone();

    // HH = final_trans(cv::Rect(cols/2, rows/2, cols/2, rows/2)).clone();
    return final_trans;
    
}

void inverse_wavelet_transfrom_1D(vector<float>& output, 
                          const vector<float>& approx,
                          const vector<float>& detail)
{
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

cv::Mat inverse_wavelet_transform_2D(cv::Mat& transformed)
{
    int rows = transformed.rows;
    int cols = transformed.cols;

    // intermediate results
    cv::Mat intermediate(rows, cols, CV_32F);
    cv::Mat image(rows, cols, CV_32F);

    // Apply inverse transform to columns
    for (int j = 0; j < cols; ++j) {
        std::vector<float> colApprox(rows / 2), colDetail(rows / 2), colOutput(rows);
        for (int i = 0; i < rows / 2; ++i) {
            colApprox[i] = transformed.at<float>(i, j);
            colDetail[i] = transformed.at<float>(i + rows / 2, j);
        }

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

        inverse_wavelet_transfrom_1D(rowOutput, rowApprox, rowDetail);

        for (int j = 0; j < cols; ++j) {
            image.at<float>(i, j) = rowOutput[j];
        }
    }

    return image;

} 

int main()
{
    std::cout << "Testing the template\n";

    cv::Mat image = cv::imread("../docs/cameraman.png");

    if (image.empty()){
        std::cerr << "Error : could not open the image" << std::endl;
    }

    image.convertTo(image, CV_32F, 1.0 / 255.0);

    cv::imshow("Cameraman image", image);

    cv::Mat LL, LH, HL, HH;
    cv::Mat WaveTransform = wavelet_transform_2D(image, 3);

    cv::imshow("After transformation", WaveTransform);

    // cv::Mat recon = inverse_wavelet_transform_2D(WaveTransform);

    // cv::imshow("Reconstructed image", recon);

    cv::waitKey(0);
    return 0;
}