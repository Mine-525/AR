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


    int subpixSampleSafe(Mat pSrc, PVector p){
        int x = (int)(floor(p.x));
        int y = (int)(floor(p.y));

        if(x < 0 || x >= pSrc.cols() - 1 || y < 0 || y >= pSrc.rows() - 1) return 127;

        int dx = (int)(256 * (p.x - floor(p.x)));
        int dy = (int)(256 * (p.y - floor(p.y)));

        int i = (int)(pSrc.get(y, x)[0]);
        int ix = (int)(pSrc.get(y, x+1)[0]);
        int iy = (int)(pSrc.get(y+1, x)[0]);
        int ixy = (int)(pSrc.get(y+1, x+1)[0]);

        int a = i + ((dx * (ix - i)) >> 8);
        int b = iy + ((dx * (ixy - iy)) >> 8);

        return a + ((dy * (b - a)) >> 8);
    }

    Point[] Pvec2Parray(PVector[] v){
        Point[] p = new Point[v.length];
        for (int i = 0; i < v.length; i++){
            p[i] = new Point(v[i].x, v[i].y);
        }
        return p;
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

            int kNumOfEdgepoints = 7;

            // circles on vertex and edge
            for (int i = 0; i < kNumOfCorners; i++){
                // vertex
                stroke(200, 0, 0);
                fill(200, 0, 0);
                circle((float)p[i].x, (float)p[i].y, kCircleSize);

                Point[] approx_points = contour_approx.toArray();
                PVector pa = OpenCV.pointToPVector(approx_points[(i+1)%kNumOfCorners]);
                PVector pb = OpenCV.pointToPVector(approx_points[i]);
                PVector kEdgeDirectionVec = PVector.div(PVector.sub(pa, pb), kNumOfEdgepoints);
                float kEdgeDirectionVecNorm = kEdgeDirectionVec.mag();

                // stripe size
                int stripe_length = (int)(0.4 * kEdgeDirectionVecNorm);
                if(stripe_length < 5) stripe_length = 5;
                stripe_length |= 1; // make odd


                int kStop = stripe_length >> 1; //floor(stripe_length/2);
                int kSrart = -kStop;
                int kStripeWidth = 3;
                Size kStripeSize = new Size(kStripeWidth, stripe_length);
                println(stripe_length);

                //Direction vectors
                PVector kStripeVecX = PVector.div(kEdgeDirectionVec, kEdgeDirectionVecNorm);
                PVector kStripeVecY = new PVector(-kStripeVecX.y, kStripeVecX.x);

                // edge, for each stripe
                for(int j = 1; j < kNumOfEdgepoints; j++){
                    PVector edge_point = PVector.add(pb, PVector.mult(kEdgeDirectionVec, j));
                    stroke(40, 90, 216);
                    strokeWeight(kStripeWidth);
                    // fill(40, 90, 216);
                    circle((float)edge_point.x, (float)edge_point.y, (float)line_width);

                    line(edge_point.x + (-stripe_length)*kStripeVecY.x,
                    edge_point.y + (-stripe_length)*kStripeVecY.y,
                    edge_point.x + (stripe_length)*kStripeVecY.x,
                    edge_point.y + (stripe_length)*kStripeVecY.y);

                    // generating stripe matrix
                    Mat stripe_image = new Mat(kStripeSize, CvType.CV_8UC1);
                    for (int m = -1; m <= 1; m++){
                        // stripe length
                        for (int n = kSrart; n <= kStop; n++){
                            PVector subpixel = PVector.add(
                                PVector.add(edge_point, PVector.mult(kStripeVecX, m)),
                                PVector.mult(kStripeVecY, n)
                            );

                            // fetch subpixel value
                            int kSubpixelValue = subpixSampleSafe(image_gray, subpixel);
                            int kStripeX = m + 1;
                            int kStripeY = n + (stripe_length >> 1);
                            stripe_image.put(kStripeY, kStripeX, kSubpixelValue);
                        }
                    }

                    // use sobel operator on stripe
                    double[] sobel_values = new double[stripe_length - 2];
                    double sobel_max = 0;
                    int idx_max = -1;

                    for (int n = 1; n < stripe_length - 1; n++){
                        byte[] pbyte = new byte[3];
                        stripe_image.get(n - 1, 0, pbyte);
                        double r1 = -(pbyte[0] & 0xFF) - 2.0 * (pbyte[1] & 0xFF) - (pbyte[2] & 0xFF);

                        stripe_image.get(n + 1, 0, pbyte);
                        double r3 = (pbyte[0] & 0xFF) + 2.0 * (pbyte[1] & 0xFF) + (pbyte[2] & 0xFF);

                        sobel_values[n-1] = -(r1 + r3);
                        if(sobel_max < sobel_values[n-1]){
                            sobel_max = sobel_values[n-1];
                            idx_max = n-1;
                        }
                    }

                    double pos, y0, y1, y2;
                    y0 = (idx_max <= 0) ? 0: sobel_values[idx_max-1]
                    y1 = sobel_values[idx_max];
                    y2 = (idx_max >= stripe-3) ? 0 : sobel_values[idx_max+1];

                    // parabolic curve
                    pos = (y2 - y0) / (4 * y1 - 2 * y0 - 2 * y2);

                    int max_idxshift = idx_max - (stripe_length >> 1);
                    PVector edge_center = PVector.add(edge_point, PVector.mult(kStripeVecY, max_idxshift + (float)pos))
                    

                }

                





                
            }
        }



        
    }
}