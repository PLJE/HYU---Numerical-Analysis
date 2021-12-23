#include <opencv2/imgcodecs.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc.hpp>
#include <string>
#include <iostream>
#include <stdio.h>
#include <math.h>

using namespace cv;
using namespace std;

int main(int ac, char** av) { //bilinear interpolation 

	Mat img = imread("dog.png");

	double M_rate = 1.0;
	double N_rate = 1.5; //MXN

	Mat resample(img.rows * M_rate, img.cols * N_rate, CV_8UC3);
	//CV_8U 는 이루어진 data type 8bit , unsigned
//  를 의미한다
	//C3 : unsigned 8bit * 3(R,G,B) 

	for (int y = 0; y < resample.rows - 1; y++) { //height
		for (int x = 0; x < resample.cols - 1; x++) { //width

			int px = (int)(x / N_rate);
			int py = (int)(y / M_rate);

			if (px >= img.cols - 1 || py >= img.rows - 1) break;

			double fx1 = (double)x / (double)N_rate - (double)px;
			double fx2 = 1 - fx1;
			double fy1 = (double)y / (double)M_rate - (double)py;
			double fy2 = 1 - fy1;

			double w1 = fx2 * fy2;
			double w2 = fx1 * fy2;
			double w3 = fx2 * fy1;
			double w4 = fx1 * fy1;

			Vec3b p1 = img.at<Vec3b>(py, px);
			Vec3b p2 = img.at<Vec3b>(py, px + 1);
			Vec3b p3 = img.at<Vec3b>(py + 1, px);
			Vec3b p4 = img.at<Vec3b>(py + 1, px + 1);

			resample.at<Vec3b>(y, x) = (w1 * p1) + (w2 * p2) + (w3 * p3) + (w4 * p4);
		}
	}

	imshow("before resampling", img);
	imshow("after resampling", resample);

	cout << "M:N ( " << M_rate << " : " << N_rate << " ) \n";
	cout << "before resampling : " << img.rows << " X " << img.cols << "\n";
	cout << "after resampling: " << resample.rows << " X " << resample.cols;

	waitKey(0);

	return 0;
}
