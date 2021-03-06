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
    PMatrix3D pose;
    Mat marker_image;
    Mat transform;
    int id;

    Marker() {
        pose = new PMatrix3D();
        marker_image = new Mat(kNumMarkerPxl, kNumMarkerPxl, CvType.CV_8UC1);
        transform = new Mat(3, 3, CvType.CV_64FC1);
    }

    void print_matrix() {
        pose.print();
    }
}

class MarkerTracker {
    boolean check_ID = false;

    int thresh;     // Threshold: gray to mono for whole image
    int bw_thresh;  // Threshold for ID image of gray marker
    double kMarkerSizeLength;
    double kMarkerSizeMin;
    double kMarkerSizeMax;

    Mat image_bgr;
  Mat image_gray;
  Mat image_gray_filtered;


    MarkerTracker(double _kMarkerSizeLength) {
        thresh = 80;
        bw_thresh = 150;
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

    PVector[] Parray2Pvector(Point[] p){
        PVector[] v = new PVector[p.length];
        for (int i = 0; i < p.length; i++){
            v[i] = new PVector((float)p[i].x, (float)p[i].y);
        }
        return v;
    }

    Point[] get_intersect(Mat[] line){
        Point[] intersections = new Point[kNumOfCorners];
        // get line parameters
        double[][] line_params = new double[kNumOfCorners][3];
        for (int i = 0; i < kNumOfCorners; i++){
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


    int[] get_ID(Marker marker){
        int[] result;
        result = new int[2];
        result[0] = (int)1e9;
        result[1] = -1;
        int sec_size = marker.marker_image.rows()/(kNumMarkRange+2);
        Mat image_refine = new Mat(kNumMarkRange, kNumMarkRange, CvType.CV_8UC1);
        int offset = 10;

        // get small marker image
        for (int i = 1; i <= kNumMarkRange; i++){
            for(int j = 1; j <= kNumMarkRange; j++){
                Mat submat = marker.marker_image.submat(sec_size * j + offset, sec_size * (j+1)-offset, sec_size * i+offset, sec_size * (i+1)-offset);
                double sum = Core.sumElems(submat).val[0];
                int ave = (int)(sum / (offset * offset));
                
                if(ave >= 120){
                    image_refine.put(i-1, j-1, 0);
                }else{
                    image_refine.put(i-1, j-1, 1);
                }
            }
        }

        // calculate ID from eacn direction and select the smallest one
        int base = (int)Math.pow(2, image_refine.rows()); // for convert id
        int direction = -1;
        
        for (int i = 0; i < 4; i++){
            int id = 0;
            for (int j = 0; j < image_refine.rows(); j++){
                int row = 0;
                for (int k = 0; k < image_refine.cols(); k++){
                    row = row + (int)image_refine.get(j, k)[0] * (int)Math.pow(2, image_refine.cols()-1-k);
                }
                id = id + row * (int)Math.pow(base, image_refine.rows()-1-j);
            }

            if(id < result[0]){
                result[0] = id;
                result[1] = i;
            }
            image_refine = rotate_mat(image_refine);
            
        }
        if(result[0] == (int)Math.pow(base, kNumMarkRange) - 1) result[0] = 0;


        return result;
    }

    Mat rotate_mat(Mat mat){
        Mat mat_rot = new Mat(mat.cols(), mat.rows(), mat.type());

        for(int i = 0; i < mat.rows(); i++){
            for(int j = 0; j < mat.cols(); j++){
                mat_rot.put(j, (mat_rot.cols() - 1) - i, mat.get(i, j));
            }
        }
        return mat_rot;
    }



  void findMarker(ArrayList<Marker> markers) {
        markers.clear();

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
        int marker_cnt = 0;



        ArrayList<MatOfPoint> contours = new ArrayList<MatOfPoint>();
        Mat hierarchy = new Mat();

        image_gray = OpenCV.imitate(opencv.getGray());
        opencv.getGray().copyTo(image_gray);
        image_gray_filtered = OpenCV.imitate(opencv.getGray());

        PImage dst = createImage(image_width, image_height, ARGB);
        opencv.toPImage(image_bgr, dst);
        
        image(dst, 0, 0);
        // beginShape();
        //     texture(dst);
        //     vertex(0,0,0,0,0);
        //     vertex(width,0,0,width,0);
        //     vertex(width,height,0,width,height);
        //     vertex(0,height,0,0,height);
        // endShape();
        

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
                vertex((float)p[i].x, (float)p[i].y, 0);
            }
            endShape(CLOSE);

            int kNumOfEdgepoints = 7;
            Mat[] line_parameters = new Mat[kNumOfCorners];
            boolean edge_detection = true;

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

                    line(edge_point.x + (-stripe_length)*kStripeVecY.x,
                    edge_point.y + (-stripe_length)*kStripeVecY.y, 0,
                    edge_point.x + (stripe_length)*kStripeVecY.x,
                    edge_point.y + (stripe_length)*kStripeVecY.y, 0);

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


                    if(idx_max == -1){
                        edge_detection = false;
                        continue;
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
                if(!edge_detection) continue;

                MatOfPoint2f mat = new MatOfPoint2f();

                // convert edge point to mat
                mat.fromArray(Pvec2Parray(edge_points));

                // fit line
                line_parameters[i] = new Mat();
                Imgproc.fitLine(mat, line_parameters[i], Imgproc.CV_DIST_L2, 0, 0.01, 0.01);

                // // display modified line
                // double vx = line_parameters[i].get(0,0)[0], vy = line_parameters[i].get(1,0)[0], x0 = line_parameters[i].get(2,0)[0], y0 = line_parameters[i].get(3,0)[0];
                // double lefty = ((-x0*vy/vx) + y0);
                // double righty = (((image_width-x0)*vy/vx)+y0);

                // line((float)(image_width-1), (float)righty, (float)0, (float)lefty);
            }

            if(!edge_detection) continue;

            // compute corners as intersections of sides
            Point[] intersections = new Point[kNumOfCorners];
            intersections = get_intersect(line_parameters);
            
            if(marker_cnt == 0){
                for (int i = 0; i < kNumOfCorners; i++){
                    switch (i) {
                        case 0:
                            stroke(200, 0, 0);
                            fill(200, 0, 0);
                            circle((float)intersections[i].x, (float)intersections[i].y, kCircleSize);
                            break;
                        case 1:
                            stroke(200, 200, 0);
                            fill(200, 200, 0);
                            circle((float)intersections[i].x, (float)intersections[i].y, kCircleSize);
                            break;
                        case 2:
                            stroke(0, 200, 0);
                            fill(0, 200, 0);
                            circle((float)intersections[i].x, (float)intersections[i].y, kCircleSize);
                            break;
                        case 3:
                            stroke(0, 200, 200);
                            fill(0, 200, 200);
                            circle((float)intersections[i].x, (float)intersections[i].y, kCircleSize);
                            break;
                    }
                }
            }
            

            // setting src image position, e.g. marker position in the image, for transform
            double src_point[] = new double[2*kNumOfCorners];
            for(int i = 0; i < kNumOfCorners; i++){
                src_point[2*i] = intersections[i].x;
                src_point[2*i+1] = intersections[i].y;
            }

            Mat src_point_mat = new Mat(4, 2, CvType.CV_32FC1);
            src_point_mat.put(0, 0, src_point);

            // setting dst image position for transform
            double dst_point[] = new double[]{0, 0, 0, kNumMarkerPxl, kNumMarkerPxl, kNumMarkerPxl, kNumMarkerPxl, 0};
            Mat dst_point_mat = new Mat(4, 2, CvType.CV_32FC1);
            dst_point_mat.put(0, 0, dst_point);

            // create marker instance for registration, and register its transform matrix and marker image
            Marker marker = new Marker(); 
            marker.transform = Imgproc.getPerspectiveTransform(src_point_mat, dst_point_mat);
            Imgproc.warpPerspective(image_gray, marker.marker_image, marker.transform, new Size((double)kNumMarkerPxl, (double)kNumMarkerPxl));
            
            Imgproc.threshold(marker.marker_image, marker.marker_image, bw_thresh, 255.0, Imgproc.THRESH_BINARY); // thresholding marker image

            // PImage dst2 = createImage(kNumMarkerPxl, kNumMarkerPxl, ARGB);
            // opencv.toPImage(marker.marker_image, dst2);
            // image(dst2, 0, 0);

            // get ID of marker
            int[] id_and_direction = new int[2];
            id_and_direction = get_ID(marker);
            marker.id = id_and_direction[0];
            if(marker.id == 0) continue;

            Point[] intersections_sub = new Point[intersections.length];

            // correct image and  for consistent marker direction
            if(id_and_direction[1] != 0){
                for (int i = 1; i <= id_and_direction[1]; i++){
                    System.arraycopy(intersections, 0, intersections_sub, 0, intersections.length);
                    marker.marker_image = rotate_mat(marker.marker_image); // rotate marker image

                    // shift intersections array index for estimateSquarePose
                    for (int j = 0; j < kNumOfCorners; j++){
                        intersections[(j+1)%kNumOfCorners] = new Point(intersections_sub[j].x, intersections_sub[j].y);
                    }
                }
            }

            
            



            
            
            

            // get pose of marker
            PVector[] intersections_v = Parray2Pvector(intersections);
            marker.pose = estimateSquarePose(intersections_v, kMarkerSize);
            



            markers.add(marker); // register marker
            marker_cnt = marker_cnt + 1;

            
            
        }


        
        
        marker_cnt = 0;
    }
}
