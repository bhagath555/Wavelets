#include <iostream>
#include <opencv2/opencv.hpp>

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