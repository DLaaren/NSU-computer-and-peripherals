/*It’s also important to note at this point that OpenCV
reads color images in BGR format, whereas most other
computer vision libraries use the RGB channel format order. 
So, when using OpenCV with other toolkits, don’t forget 
to swap the blue and red color channels, 
as you switch from one library to another.*/

#include <iostream>
#include <opencv2/opencv.hpp> 

using namespace std;
using namespace cv;

int main() {
    //reading an image
    Mat img_color = imread("rose.jpg", 1); //1 = default
    Mat img_gray = imread("rose.jpg", 0);
    Mat img_unchanged = imread("rose.jpg", -1);

    //displaying an image
    //1) create a window
    namedWindow("color", WINDOW_AUTOSIZE);
    namedWindow("gray", WINDOW_AUTOSIZE);
    namedWindow("unchanged", WINDOW_AUTOSIZE);
    //2) show the image inside the window
    imshow("color", img_color);
    imshow("gray", img_gray);
    imshow("unchanged", img_unchanged);
    //wait for a keystroke
    waitKey(0);
    //destroy all the windows created
    destroyAllWindows();
    return 0;
}
