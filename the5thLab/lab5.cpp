#include <iostream>
#include <ctime>
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
    int fps = 15;
    VideoWriter output("output.mp4", VideoWriter::fourcc('M', 'P', '4', 'V'),fps, frame_size);
    
    double totalTime = 0, inputTime = 0, procTime = 0, outputTime = 0;
    int frameCount = 0;

    cout << "Video capturing has been started" << endl;
    //now we need to write the video per frame
    while (web_cam.isOpened()) {
        Mat frame;
        bool isSuccess = web_cam.read(frame);
        if (!isSuccess) {
            cout << "stream disconnected" << endl;
            break;
        }
        //reading the video
        clock_t start = clock();
        web_cam >> frame;
        clock_t end = clock();
        double write_time = 1000.0 * (end - start) / CLOCKS_PER_SEC;
        totalTime += write_time;
        inputTime += write_time;

        //proccessing the video
        start = clock();
        applyColorMap(frame, frame, COLORMAP_TWILIGHT); //look colormaps in the docs
        end = clock();
        double proc_time = 1000.0 * (end - start) / CLOCKS_PER_SEC;
        totalTime += proc_time;
        procTime += proc_time;

        //writing the video
        output.write(frame);

        //showing the video
        start = clock();
        imshow("Result", frame);
        end = clock();
        double output_time = 1000.0 * (end - start) / CLOCKS_PER_SEC;
        totalTime += output_time;
        outputTime += output_time;
        frameCount++;

        int key = waitKey(20);
        if (key == 'q') {
            cout << "Key q key is pressed. Stopping the video" << endl;
            break;
        }
    }
    double procent = totalTime / 100.0;
    cout << "CPU time for reading frames: " << inputTime / procent << "%" << endl;
    cout << "CPU time for processing frames: " << procTime / procent << "%" << endl;
    cout << "CPU time for output frames: " << outputTime / procent << "%" << endl;
    cout << "Total framerate: " << ((frameCount) / (totalTime / 1000.0)) << endl;
    return 0;
}
