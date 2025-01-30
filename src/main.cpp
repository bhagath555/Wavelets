#include "wavelet_transform.hpp"

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
void wavelet_transform_1D(const vector<float>& input,
                          vector<float>& approx,
                          vector<float>& detail,
                          const vector<float>& lowpass,
                          const vector<float>& highpass) {
    if (input.empty()) {
        throw std::invalid_argument("Input vector is empty");
    }
    size_t filter_size = lowpass.size();
    size_t size = input.size() / 2;

    approx.resize(size);
    detail.resize(size);

    // Apply lowpass and highpass filters
    for (size_t i = 0; i < size; ++i) {
        approx[i] = 0;
        detail[i] = 0;

        for (size_t k = 0; k < filter_size; ++k) {
            size_t index = (2 * i + k) % input.size(); // Handle index out of range.
            approx[i] += input[index] * lowpass[k];
            detail[i] += input[index] * highpass[k];
        }
    }
}

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
cv::Mat wavelet_transform_2D(const cv::Mat& image, const int& level,
                             const vector<float>& lowpass, const vector<float>& highpass) {
    int rows = image.rows;
    int cols = image.cols;
    cv::Mat final_trans(rows, cols, CV_32F);

    if (level < 1) {
        return image;
    }

    int l = 1;
    while (l <= level) {
        if (l > 1) {
            rows /= 2;
            cols /= 2;
        }

        cv::Mat row_trans(rows, cols, CV_32F);
        cv::Mat row;

        // Apply the transformation for rows
        for (int i = 0; i < rows; ++i) {
            if (l == 1) {
                row = image.row(i);
            } else {
                row = final_trans.row(i).colRange(0, cols);
            }
            vector<float> rowVector(row.begin<float>(), row.end<float>());

            // Perform 1D transformation
            vector<float> approx, detail;
            wavelet_transform_1D(rowVector, approx, detail, lowpass, highpass);

            for (int j = 0; j < approx.size(); ++j) {
                row_trans.at<float>(i, j) = approx[j];
                row_trans.at<float>(i, j + cols / 2) = detail[j];
            }
        }

        // Apply the transformation for cols
        for (int j = 0; j < cols; ++j) {
            cv::Mat col = row_trans.col(j);
            vector<float> colVector(col.begin<float>(), col.end<float>());

            // Perform 1D transformation
            vector<float> approx, detail;
            wavelet_transform_1D(colVector, approx, detail, lowpass, highpass);

            for (int i = 0; i < approx.size(); ++i) {
                final_trans.at<float>(i, j) = approx[i];
                final_trans.at<float>(i + rows / 2, j) = detail[i];
            }
        }
        l++;
    }
    return final_trans;
}

int main() {
    // Example filters
    vector<float> haar_lowpass = {0.7071, 0.7071};
    vector<float> haar_highpass = {-0.7071, 0.7071};

    vector<float> db2_lowpass = {-0.1294, 0.2241, 0.8365, 0.4830};
    vector<float> db2_highpass = {-0.4830, 0.8365, -0.2241, -0.1294};

    vector<float> coif1_lowpass = {-0.0157, -0.0727, 0.3849, 0.8526, 0.3379, -0.0727};
    vector<float> coif1_highpass = {0.0727, 0.3379, -0.8526, 0.3849, 0.0727, -0.0157};

    vector<float> sym3_lowpass = {0.0352, -0.0854, -0.1350, 0.4599, 0.8069, 0.3327};
    vector<float> sym3_highpass = {-0.3327, 0.8069, -0.4599, -0.1350, 0.0854, 0.0352};


    // Load image
    cv::Mat image = cv::imread("../docs/cameraman.png", cv::IMREAD_GRAYSCALE);
    if (image.empty()) {
        std::cerr << "Error: Could not open the image." << std::endl;
        return -1;
    }

    if (image.rows % 2 != 0 || image.cols % 2 != 0) {
        throw std::runtime_error("Image dimensions must be even.");
    }
    image.convertTo(image, CV_32F, 1.0 / 255.0);

    cv::imshow("Original Image", image);

    // Perform wavelet transformation with Haar filters
    int level = 3;
    cv::Mat wavelet_haar = wavelet_transform_2D(image, level, haar_lowpass, haar_highpass);
    cv::imshow("Haar Wavelet Transform", wavelet_haar);

    // Perform wavelet transformation with db2 filters
    cv::Mat wavelet_db2 = wavelet_transform_2D(image, level, db2_lowpass, db2_highpass);
    cv::imshow("DB2 Wavelet Transform", wavelet_db2);

    // Perform wavelet transformation with coif1 filters
    cv::Mat wavelet_coif1 = wavelet_transform_2D(image, level, coif1_lowpass, coif1_highpass);
    cv::imshow("Coif1 velet Transform", wavelet_coif1);

    // Perform wavelet transformation with coif1 filters
    cv::Mat wavelet_sym3 = wavelet_transform_2D(image, level, sym3_lowpass, sym3_highpass);
    cv::imshow("Sym3 Wavelet Transform", wavelet_sym3);

    cv::waitKey(0);
    return 0;
}