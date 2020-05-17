import processing.video.*;

import gab.opencv.*;
import processing.video.*;

final boolean MARKER_TRACKER_DEBUG = true;

final boolean USE_SAMPLE_IMAGE = true;

// We've found that some Windows build-in cameras (e.g. Microsoft Surface)
// cannot work with processing.video.Capture.*.
// Instead we use DirectShow Library to launch these cameras.
final boolean USE_DIRECTSHOW = false;

final double kMarkerSize = 0.03; // [m]
final int thresh = 70;
final int bw_thresh = 70;
final int kNumOfCorners = 4;

Capture cam;
OpenCV opencv;

ArrayList<Marker> markers;
MarkerTracker markerTracker;

PImage img;

public void settings(){
  if(USE_SAMPLE_IMAGE){
      size(1000, 730);
  }else{
    size(640, 480);
  }
}


void setup() {
  
  if(USE_SAMPLE_IMAGE){
    img = loadImage("./marker_test.jpg");
    opencv = new OpenCV(this, "./marker_test.jpg");
    image(img, 0, 0);
  }else{
    String[] cameras = Capture.list();

    
    if (cameras.length == 0) {
      println("There are no cameras available for capture.");
      exit();
    } else {
      println("Available cameras:");
      for (int i = 0; i < cameras.length; i++) {
        println(cameras[i]);
      }
      
      // The camera can be initialized directly using an 
      // element from the array returned by list():
      cam = new Capture(this, cameras[0]);
        

      opencv = new OpenCV(this, cam.width, cam.height);   
    }
    if(!USE_DIRECTSHOW){
      cam.start();
    }
  }

  smooth();
  markerTracker = new MarkerTracker(kMarkerSize, thresh, bw_thresh);
  
}

void draw() {
  if(!USE_SAMPLE_IMAGE){
    if(cam.width <= 0 || cam.height <= 0){
      println("Incorrect capture data");
      return;
    }
    opencv.loadImage(cam);
  }

  markerTracker.findMarker(markers);

  System.gc();

}

void captureEvent(Capture c) {
  if (!USE_DIRECTSHOW && c.available())
      c.read();
}
