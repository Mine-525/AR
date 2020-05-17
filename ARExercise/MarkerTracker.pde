import gab.opencv.*;
import org.opencv.imgproc.Imgproc;

import org.opencv.core.Core;

import org.opencv.core.Mat;
import org.opencv.core.MatOfPoint;
import org.opencv.core.MatOfPoint2f;
import org.opencv.core.MatOfPoint2f;
import org.opencv.core.CvType;

import org.opencv.core.Point;
import org.opencv.core.Size;
import org.opencv.core.Rect;


class Marker {
    int code;
    float[] pose;

    Marker() {
        pose = new float[16];
    }

    void print_matrix() {
        int kSize = 4;
        for (int r = 0; r < kSize; r++) {
            for (int c = 0; c < kSize; c++) {
                println(pose[r + kSize * (c - 1)]);
            }
        }
    }
}

class MarkerTracker {
    int thresh;     // Threshold: gray to mono
    int bw_thresh;  // Threshold for gray marker to ID image
    double kMarkerSizeLength;
    double kMarkerSizeMin;
    double kMarkerSizeMax;

    Mat image_bgr;
	Mat image_gray;
	Mat image_gray_filtered;

    MarkerTracker(double _kMarkerSizeLength) {
        thresh = 80;
        bw_thresh = 100;
        kMarkerSizeLength = _kMarkerSizeLength;
        init();
    }

    MarkerTracker(double _kMarkerSizeLength, int _thresh, int _bw_thresh) {
        thresh = _thresh;
        bw_thresh = _bw_thresh;
        kMarkerSizeLength = _kMarkerSizeLength;
        init();
    }

    void init() {
        println("Startup");
    }

    void cleanup() {
        println("Finished");
    }

	void findMarker(ArrayList<Marker> markers) {
        boolean isFirstStripe = true;
        boolean isFirstMarker = true;

        image_bgr = OpenCV.imitate(opencv.getColor());
        opencv.getColor().copyTo(image_bgr);
        int image_height = image_bgr.rows();
        int image_width  = image_bgr.cols();
        int image_area = image_height * image_width;

        kMarkerSizeMax = image_area / 2;
        kMarkerSizeMin = image_area / 100;

        int kCircleSize = 10;
        int line_width = 4;



        ArrayList<MatOfPoint> contours = new ArrayList<MatOfPoint>();
        Mat hierarchy = new Mat();

        image_gray = OpenCV.imitate(opencv.getGray());
        opencv.getGray().copyTo(image_gray);
        image_gray_filtered = OpenCV.imitate(opencv.getGray());

        PImage dst = createImage(image_width, image_height, ARGB);
        opencv.toPImage(image_gray, dst);
        image(dst, 0, 0);
        

        // Thresholding and find contour
        Imgproc.threshold(image_gray, image_gray_filtered, thresh, 255.0, Imgproc.THRESH_BINARY);
        Imgproc.findContours(image_gray_filtered, contours, hierarchy, Imgproc.RETR_LIST, Imgproc.CHAIN_APPROX_SIMPLE);

        // Marker detection
        for (MatOfPoint contour : contours){
            MatOfPoint2f contour_approx = new MatOfPoint2f();
            double kEpsilon = 0.05 * Imgproc.arcLength(new MatOfPoint2f(contour.toArray()), true); // approximation accuracy

            Imgproc.approxPolyDP(new MatOfPoint2f(contour.toArray()), contour_approx, kEpsilon, true);

            // Check size
            Rect bounding_rectangle = Imgproc.boundingRect(new MatOfPoint(contour_approx.toArray()));
            double marker_size = bounding_rectangle.area();
            boolean is_contour_valid = (marker_size > kMarkerSizeMin)
                && (marker_size < kMarkerSizeMax)
                && (contour_approx.size().height == kNumOfCorners)
                && Imgproc.isContourConvex(new MatOfPoint(contour_approx.toArray()));

            // println(is_contour_valid);
            if(!is_contour_valid) continue;

            // drow shapes
            stroke(0, 130, 0);
            strokeWeight(line_width);
            noFill();
            
            // rectangle
            Point[] p = contour_approx.toArray();
            beginShape();
            for (int i = 0; i < p.length; i++){
                vertex((float)p[i].x, (float)p[i].y);
            }
            endShape(CLOSE);

            // circles on vertex and edge
            for (int i = 0; i < 4; i++){
                // vertex
                stroke(200, 0, 0);
                fill(200, 0, 0);
                circle((float)p[i].x, (float)p[i].y, kCircleSize);

                // edge
                for(int j = 1; j < 7; j++){
                    Point edge_point = new Point();
                    edge_point.x = p[i%4].x*j/7 + p[(i+1)%4].x*(7-j)/7;
                    edge_point.y = p[i%4].y*j/7 + p[(i+1)%4].y*(7-j)/7;
                    stroke(40, 90, 216);
                    fill(40, 90, 216);
                    circle((float)edge_point.x, (float)edge_point.y, (float)line_width);
                }


            }
        }



        
    }
}