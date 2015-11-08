#include "homography.h"
#include "matrix.h"

using namespace std;


void applyHomography(const Image &source, const Matrix &H, Image &out, bool bilinear) {
    // // --------- HANDOUT  PS06 ------------------------------
    // Transform image source using the homography H, and composite in onto out.
    // if bilinear == true, using bilinear interpolation. Use nearest neighbor
    // otherwise.

    Matrix inv_H = H.inverse();

    for (int a = 0; a < out.width(); a++){
        for (int b = 0; b < out.height(); b++){
            for (int c = 0; c < out.channels(); c++){
                
                Matrix homographic_points = Matrix(3,1);
                homographic_points(0,0) = a;
                homographic_points(1,0) = b;
                homographic_points(2,0) = 1;

                Matrix soln = inv_H * homographic_points;

                float new_x = soln(0,0)/soln(2,0);
                float new_y = soln(1,0)/soln(2,0);

                if (new_x >= 0 && new_y >= 0 && new_x < source.width() - 1&& new_y < source.height() - 1){
                    float pixel_value;
                    if (bilinear){
                        pixel_value = interpolateLin(source, new_x, new_y, c, false);
                        //interpolateLin(const Image &im, float x, float y, int z, bool clamp)
                    }
                    else {
                        pixel_value = source(round(new_x), round(new_y), c);
                    }

                    out(a, b, c) = pixel_value;

                }

            }
        }
    }
}




Matrix computeHomography(const CorrespondencePair correspondences[4]) {
    // --------- HANDOUT  PS06 ------------------------------
    // Compute a homography from 4 point correspondences.
    Matrix A = Matrix::Zero(9,9);

    //cout << "hello!" << endl;

    for (int i = 0; i < 4; i++){
        Vec3f point1 = correspondences[i].point1;
        Vec3f point2 = correspondences[i].point2;
        
        //std::cout << "adding to rows: " << 2*i << " and " << 2*i + 1 << endl;

        float x = point1.x();
        float y = point1.y();
        float x_prime = point2.x();
        float y_prime = point2.y();

        A(2*i, 0) = y;
        A(2*i, 1) = x;
        A(2*i, 2) = 1;
        A(2*i, 3) = 0;
        A(2*i, 4) = 0;
        A(2*i, 5) = 0;
        A(2*i, 6) = -y*y_prime;
        A(2*i, 7) = -x*y_prime;
        A(2*i, 8) = -y_prime;

        A(2*i + 1, 0) = 0;
        A(2*i + 1, 1) = 0;
        A(2*i + 1, 2) = 0;
        A(2*i + 1, 3) = y;
        A(2*i + 1, 4) = x;
        A(2*i + 1, 5) = 1;
        A(2*i + 1, 6) = -x_prime*y;
        A(2*i + 1, 7) = -x_prime*x;
        A(2*i + 1, 8) = -x_prime;

    }

    //...........I think these should all be ones?

    A(8,8) = 1;
    
    Matrix b = Matrix::Zero(9,1);
    b(8,0) = 1;

    Matrix x = A.fullPivLu().solve(b);

    //std::cout << A << endl;
    //cout << "x: " << endl << x << endl;

    Matrix H = Matrix::Zero(3,3);
    H(0,0) = x(0,0);
    H(0,1) = x(1,0);
    H(0,2) = x(2,0);
    H(1,0) = x(3,0);
    H(1,1) = x(4,0);
    H(1,2) = x(5,0);
    H(2,0) = x(6,0);
    H(2,1) = x(7,0);
    H(2,2) = x(8,0);

    Matrix H_prime = Matrix::Zero(3,3);
    H_prime(0,0) = H(1,1);
    H_prime(0,1) = H(1,0);
    H_prime(0,2) = H(1,2);
    H_prime(1,0) = H(0,1);
    H_prime(1,1) = H(0,0);
    H_prime(1,2) = H(0,2);
    H_prime(2,0) = H(2,1);
    H_prime(2,1) = H(2,0);
    H_prime(2,2) = H(2,2);    
    
    //cout << "H is: " << endl << H << endl;

    return H_prime;
}


BoundingBox computeTransformedBBox(int imwidth, int imheight, Matrix H) {
    // --------- HANDOUT  PS06 ------------------------------
    // Predict the bounding boxes that encompasses all the transformed
    // coordinates for pixels frow and Image with size (imwidth, imheight)
    Vec3f lower_left(0, imheight, 1);
    Vec3f lower_right(imwidth, imheight, 1);
    Vec3f upper_left(0,0,1);
    Vec3f upper_right(imwidth, 0, 1);

    Vec3f transformed_lower_left = H*lower_left;
    Vec3f transformed_lower_right = H*lower_right;
    Vec3f transformed_upper_left = H*upper_left;
    Vec3f transformed_upper_right = H*upper_right;

    float lower_left_x = transformed_lower_left.x()/transformed_lower_left.z();
    float lower_right_x = transformed_lower_right.x()/transformed_lower_right.z();
    float upper_left_x = transformed_upper_left.x()/transformed_upper_left.z();
    float upper_right_x = transformed_upper_right.x()/transformed_upper_right.z();

    float lower_left_y = transformed_lower_left.y()/transformed_lower_left.z();
    float lower_right_y = transformed_lower_right.y()/transformed_lower_right.z();
    float upper_left_y = transformed_upper_left.y()/transformed_upper_left.z();
    float upper_right_y = transformed_upper_right.y()/transformed_upper_right.z();

    return BoundingBox(min(lower_left_x, upper_left_x), max(lower_right_x, upper_right_x), min(upper_left_y, upper_right_y), max(lower_left_y, lower_right_y));
}


BoundingBox bboxUnion(BoundingBox B1, BoundingBox B2) {
    // --------- HANDOUT  PS06 ------------------------------
    // Compute the bounding box that tightly bounds the union of B1 an B2.
    return BoundingBox(min(B1.x1, B2.x1), max(B1.x2,B2.x2), min(B1.y1, B2.y1), max(B1.y2, B2.y2));
    
}


Matrix makeTranslation(BoundingBox B) {
    // --------- HANDOUT  PS06 ------------------------------
    // Compute a translation matrix (as a homography matrix) that translates the
    // top-left corner of B to (0,0).
    float x_translate = -B.x1;
    float y_translate = -B.y1;
    Matrix translation_matrix = Matrix::Identity(3,3);
    translation_matrix(0,2) = x_translate;
    translation_matrix(1,2) = y_translate;
    return translation_matrix;
}


Image stitch(const Image &im1, const Image &im2, const CorrespondencePair correspondences[4]) {
    // --------- HANDOUT  PS06 ------------------------------
    // Transform im1 to align with im2 according to the set of correspondences.
    // make sure the union of the bounding boxes for im2 and transformed_im1 is
    // translated properly (use makeTranslation)
    Matrix H = computeHomography(correspondences);
    BoundingBox B_1 = computeTransformedBBox(im1.width(), im1.height(), H);
    BoundingBox B_2 = computeTransformedBBox(im2.width(), im2.height(), Matrix::Identity(3,3));
    BoundingBox B = bboxUnion(B_1,B_2);
    Matrix T = makeTranslation(B);

    std::cout << "HELO HUNAM" << endl;
    Image out(B.x2 - B.x1, B.y2 - B.y1, im1.channels());    
    std::cout << "HELO HUNAM 1" << endl;

    applyHomographyFast(im2, T, out, true);
    std::cout << "HELO HUNAM 4" << endl;

    applyHomographyFast(im1, T*H, out, true);
    std::cout << "HELO HUNAM 2" << endl;

    return out;
}

// debug-useful
Image drawBoundingBox(const Image &im, BoundingBox bbox) {
    // // --------- HANDOUT  PS06 ------------------------------
    //  ________________________________________
    // / Draw me a bounding box!                \
    // |                                        |
    // | "I jumped to my                        |
    // | feet, completely thunderstruck. I      |
    // | blinked my eyes hard. I looked         |
    // | carefully all around me. And I saw a   |
    // | most extraordinary small person, who   |
    // | stood there examining me with great    |
    // | seriousness."                          |
    // \              Antoine de Saint-Exupery  /
    //  ----------------------------------------
    //         \   ^__^
    //          \  (oo)\_______
    //             (__)\       )\/\
    //                 ||----w |
    //                 ||     ||

    Image output = im;
    std::cout << "x1 is: " << bbox.x1 << " x2: " << bbox.x2 << " y1: " << bbox.y1 << " y2: " << bbox.y2 << endl;
    for (int i = bbox.x1; i < bbox.x2; i++){
        output(i,bbox.y1) = 1;
        output(i,bbox.y2) = 1;
    }


    for (int j = bbox.y1; j < bbox.y2; j++){
        output(bbox.x1,j) = 1;
        output(bbox.x2,j) = 1;
    }

    return output;
}

void applyHomographyFast(const Image &source, const Matrix &H, Image &out, bool bilinear) {
    // // --------- HANDOUT  PS06 ------------------------------
    // Same as apply but change only the pixels of out that are within the
    // predicted bounding box (when H maps source to its new position).
    BoundingBox B = computeTransformedBBox(source.width(), source.height(), H);

    std::cout << "In applyHomographyFast " << endl;

    Matrix inv_H = H.inverse();
    for (int a = B.x1; a < B.x2 + 1; a++){
        for (int b = B.y1; b < B.y2 + 1; b++){
            for (int c = 0; c < out.channels(); c++){
                
                Matrix homographic_points = Matrix(3,1);
                homographic_points(0,0) = a;
                homographic_points(1,0) = b;
                homographic_points(2,0) = 1;

                Matrix soln = inv_H * homographic_points;

                float new_x = soln(0,0)/soln(2,0);
                float new_y = soln(1,0)/soln(2,0);

                if (new_x >= 0 && new_y >= 0 && new_x < source.width() - 1 && new_y < source.height() - 1){
                    float pixel_value;
                    if (bilinear){
                        pixel_value = interpolateLin(source, new_x, new_y, c, false);
                        //interpolateLin(const Image &im, float x, float y, int z, bool clamp)
                    }
                    else {
                        pixel_value = source(round(new_x), round(new_y), c);
                    }

                    out(a, b, c) = pixel_value;

                }

            }
        }
    }

}
