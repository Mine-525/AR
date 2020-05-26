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
    Mat mark_image;
    Mat transform;

    Marker() {
        pose = new float[16];
        mark_image = new Mat(kNumMarkerPxl, kNumMarkerPxl, CvType.CV_8UC1);
        transform = new Mat(3, 3, CvType.CV_64FC1);
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
    boolean check_ID = false;

    int thresh;     // Threshold: gray to mono
    int bw_thresh;  // Threshold for gray marker to ID image
    double kMarkerSizeLength;
    double kMarkerSizeMin;
    double kMarkerSizeMax;

    Mat image_bgr;
	Mat image_gray;
	Mat image_gray_filtered;

    ArrayList<Marker> marker_list = new ArrayList<Marker>(); // list for detected markers

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

    Point[] get_intersect(Mat[] line){
        Point[] intersections = new Point[kNumOfCorners];
        // get line parameters
        double[][] line_params = new double[kNumOfCorners][3];
        for (int i = 0; i < kNumOfCorners; i++){
            // System.out.println(line[i].dump());
            double vx = line[i].get(0,0)[0], vy = line[i].get(1,0)[0], x0 = line[i].get(2,0)[0], y0 = line[i].get(3,0)[0];
            line_params[i][0] = -vy;
            line_params[i][1] = vx;
            line_params[i][2] = vx * y0 - vy * x0;

        }

        // compute intersects
        for (int i = 0; i < kNumOfCorners; i++){
            double a = line_params[i][0], b = line_params[i][1], e = line_params[i][2];
            double c = line_params[(i+1)%4][0], d = line_params[(i+1)%4][1], f = line_params[(i+1)%4][2];

            double x = (d * e - b * f) / (a * d - b * c), y = (a * f - e * c) / (a * d - b * c);
            intersections[i] = new Point(x, y);
        }

        return intersections;
    }

    void showID(){
        if(check_ID){
            PImage mini_ID = createImage(200, 200, ARGB);
            Marker marker_1st = marker_list.get(0);
            opencv.toPImage(marker_1st.mark_image, mini_ID);
            System.out.println(marker_1st.mark_image.dump());
            image(mini_ID, 0, 0);
        }
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
            Mat[] line_parameters = new Mat[kNumOfCorners];

            // circles on vertex and edge
            for (int i = 0; i < kNumOfCorners; i++){
                // vertex
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

                //Direction vectors
                PVector kStripeVecX = PVector.div(kEdgeDirectionVec, kEdgeDirectionVecNorm);
                PVector kStripeVecY = new PVector(-kStripeVecX.y, kStripeVecX.x);

                PVector[] edge_points = new PVector[kNumOfEdgepoints-1];
                Mat[] lines = new Mat[kNumOfCorners];

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
                    y0 = (idx_max <= 0) ? 0: sobel_values[idx_max-1];
                    y1 = sobel_values[idx_max];
                    y2 = (idx_max >= stripe_length-3) ? 0 : sobel_values[idx_max+1];

                    // parabolic curve
                    pos = (y2 - y0) / (4 * y1 - 2 * y0 - 2 * y2);

                    int max_idxshift = idx_max - (stripe_length >> 1);
                    
                    edge_points[j-1] = PVector.add(edge_point, PVector.mult(kStripeVecY, max_idxshift + (float)pos));

                }

                MatOfPoint2f mat = new MatOfPoint2f();

                // convert edge point to mat
                mat.fromArray(Pvec2Parray(edge_points));

                // fit line
                line_parameters[i] = new Mat();
                Imgproc.fitLine(mat, line_parameters[i], Imgproc.CV_DIST_L2, 0, 0.01, 0.01);

                // double vx = line_parameters[i].get(0,0)[0], vy = line_parameters[i].get(1,0)[0], x0 = line_parameters[i].get(2,0)[0], y0 = line_parameters[i].get(3,0)[0];
                // double lefty = ((-x0*vy/vx) + y0);
                // double righty = (((image_width-x0)*vy/vx)+y0);

                // line((float)(image_width-1), (float)righty, (float)0, (float)lefty);
            }

            // compute corners as intersections of sides
            Point[] intersections = new Point[kNumOfCorners];
            intersections = get_intersect(line_parameters);
            
            stroke(200, 0, 0);
            fill(200, 0, 0);
            for (int i = 0; i < kNumOfCorners; i++){
                circle((float)intersections[i].x, (float)intersections[i].y, kCircleSize);
            }

            // setting src image position for transform
            double src_point[] = new double[2*kNumOfCorners];
            for(int i = 0; i < kNumOfCorners; i++){
                src_point[2*i] = intersections[i].x;
                src_point[2*i+1] = intersections[i].y;
            }

            Mat src_point_mat = new Mat(4, 2, CvType.CV_32FC1);
            src_point_mat.put(0, 0, src_point);
            // println(src_point_mat.dump());

            // setting dst image position for transform
            double dst_point[] = new double[]{0, 0, 0, 200, 200, 200, 200, 0};
            Mat dst_point_mat = new Mat(4, 2, CvType.CV_32FC1);
            dst_point_mat.put(0, 0, dst_point);
            // println(dst_point_mat.dump());


            Marker marker = new Marker();
            marker.transform = Imgproc.getPerspectiveTransform(src_point_mat, dst_point_mat);
            // println(marker.transform.dump());
            Imgproc.warpPerspective(image_gray, marker.mark_image, marker.transform, new Size((double)kNumMarkerPxl, (double)kNumMarkerPxl));
            marker_list.add(marker);
            if(!check_ID){
                check_ID = true;
            }
            Imgproc.threshold(marker.mark_image, marker.mark_image, bw_thresh, 255.0, Imgproc.THRESH_BINARY);
            
        }



        
    }
}