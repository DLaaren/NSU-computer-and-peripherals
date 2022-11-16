#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

int main() {
    //reading a file
    VideoCapture video("street.mp4");
    //check if the file opened
    if (!video.isOpened()) {
        cout << "not opend" << endl;
        return 0;
    }
    //obtain fps and frame count by get() method
    int fps = video.get(CAP_PROP_FPS);
    cout << "Frame Rate : " << fps << endl;

    int frame_count = video.get(CAP_PROP_FRAME_COUNT);
    cout << "Frame Count : " << frame_count << endl;

    while (video.isOpened()) {
        //initialize fram matrix 
        Mat frame;
        bool isSuccess = video.read(frame);
        if (isSuccess) {
            imshow("Frame", frame);
        }
        else {
            cout << "video ends" << endl;
        }
        int key = waitKey(20);
        if (key == 'q') {
            cout << "q key is pressed. Stopping the video" << endl;
            break;
        }
    }
    return 0;
}