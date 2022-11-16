#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

int main() {
    //index for the camera is 0
    //for more than one camera -- each camera incremented (e.g. 1,2,etc)
    VideoCapture web_cam(0);

    //obtain frame size
    int frame_width = static_cast<int>(web_cam.get(CAP_PROP_FRAME_WIDTH));
    int frame_height = static_cast<int>(web_cam.get(CAP_PROP_FRAME_HEIGHT));
    Size frame_size(frame_width, frame_height);
    int fps = 30;
    VideoWriter output("output.mp4", VideoWriter::fourcc('M', 'P', '4', 'V'),fps, frame_size);
    
    //now we need to write the video per frame
    while (web_cam.isOpened()) {
        Mat frame;
        bool isSuccess = web_cam.read(frame);
        if (!isSuccess) {
            cout << "stream disconnected" << endl;
            break;
        }
        else {
            output.write(frame);
            imshow("Frame", frame);
            int key = waitKey(20);
            if (key == 'q') {
                cout << "Key q key is pressed by the user. Stopping the video" << endl;
                break;
            }
        }
    }
    web_cam.release();
    output.release();
    return 0;
}
