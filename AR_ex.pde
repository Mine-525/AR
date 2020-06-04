import processing.video.*;

import gab.opencv.*;
import processing.video.*;

final boolean MARKER_TRACKER_DEBUG = true;

final boolean USE_SAMPLE_IMAGE = true;

// We've found that some Windows build-in cameras (e.g. Microsoft Surface)
// cannot work with processing.video.Capture.*.
// Instead we use DirectShow Library to launch these cameras.
final boolean USE_DIRECTSHOW = false;

final float kMarkerSize = 0.03; // [m]
final int thresh = 70;
final int bw_thresh = 170;
final int kNumOfCorners = 4;
final int kNumMarkerPxl = 210;
final int kNumMarkRange = 4;
final float fov = 45;
final float light_value = 120;


Capture cam;
OpenCV opencv;


MarkerTracker markerTracker;

PImage img;

public void settings(){
  if(USE_SAMPLE_IMAGE){
      size(1280, 720, P3D);
  }else{
    size(640, 480, P3D);
  }
}


void setup() {
  
  if(USE_SAMPLE_IMAGE){
    img = loadImage("./marker_test2.jpg");
    opencv = new OpenCV(this, "./marker_test2.jpg");
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

  // Get viewing conversion matrix M_v 
  PMatrix3D cam_mat = ((PGraphicsOpenGL)g).camera; 
 
  // set M_v as identity 
  cam_mat.reset(); 
  frameRate(60);


}

void draw() {
  ArrayList<Marker> markers = new ArrayList<Marker>();
  if(!USE_SAMPLE_IMAGE){
    if(cam.width <= 0 || cam.height <= 0){
      println("Incorrect capture data");
      return;
    }
    opencv.loadImage(cam);
  }

  // use orthographic camera 
  ortho(float(-width/2), float(width/2), float(-height/2), float(height/2), 0.0, (height/2) / tan(radians(fov))); 
  pushMatrix();
    println("image : " + modelZ(0, 0, 0));
    translate(float(-width/2), float(-height/2), -(height/2) / tan(radians(fov))); // z : -360
    println("image translate: " + modelZ(0, 0, 0));
    markerTracker.findMarker(markers);
  popMatrix();

  // use perspective camera 
  perspective(radians(fov), float(width) / float(height), 0.01, 1000.0); 

  pushMatrix();
    // set snowman's position
    println("snowman : " + modelZ(0, 0, 0));
    
    translate(0, 2, -10); 
    println("snowman translate : " + modelZ(0, 0, 0));
    directionalLight(light_value, light_value, light_value, -10, 10, 0); 
     
    rotateY(radians(frameCount * 2)); 
    rotateX(radians(0)); 
    rotateZ(radians(0)); 
 
    // draw snowman 
    snowman(1);
  popMatrix();

  noLights();

  System.gc();

}

void captureEvent(Capture c) {
  if (!USE_DIRECTSHOW && c.available())
      c.read();
}
