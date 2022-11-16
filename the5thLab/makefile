CC=g++
CFLAGS=-std=c++20
INCLUDEDIR=/usr/local/include
LIBDIR=/usr/local/lib
LIBS=-lopencv_calib3d -lopencv_core -lopencv_dnn -lopencv_features2d -lopencv_flann -lopencv_gapi -lopencv_highgui -lopencv_imgcodecs -lopencv_imgproc -lopencv_ml -lopencv_objdetect -lopencv_photo -lopencv_stitching -lopencv_video -lopencv_videoio
SRC=lab5.cpp
DEST=lab


all:
	$(CC) $(CFLAGS) -I$(INCLUDEDIR) -L$(LIBDIR) $(LIBS) $(SRC) -o $(DEST) 