#include "blending.h"
#include "matrix.h"
#include <ctime>

using namespace std;

Image blendingweight(int imwidth, int imheight) {
    // // --------- HANDOUT  PS07 ------------------------------
    Image output(imwidth, imheight, 1);

    for (int i = 0; i < imwidth; i++){
        for (int j = 0; j < imheight; j++){
            output(i,j) =  (1 - fabs(imwidth/2  - i)/(imwidth/2)) * ( 1 - fabs(imheight/2  - j)/(imheight/2));
        }
    }
    
    return output;
}

//  ****************************************************************************
//  * blending related functions re-written from previous asasignments
//  ****************************************************************************

// instead of writing source in out, *add* the source to out based on the weight
// so out(x,y) = out(x, y) + weight * image
void applyhomographyBlend(const Image &source, const Image &weight, Image &out, Matrix &H, bool bilinear) {
    // // --------- HANDOUT  PS07 ------------------------------

    Matrix inv_H = H.inverse();

    cout << "TODO: WEIRD IMAGE BOUNDS THING " << endl;

    cout << "weight image width: " << weight.width() << " height: " << weight.height() << endl;

    cout << "out image width: " << out.width() << " height: " << out.height() << endl;


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

                if (new_x >= 0 && new_y >= 0 && new_x < source.width() - 1 && new_y < source.height() - 1){
                    float pixel_value;
                    if (bilinear){
                        pixel_value = interpolateLin(source, new_x, new_y, c, false);
                        //interpolateLin(const Image &im, float x, float y, int z, bool clamp)
                    }
                    else {
                        pixel_value = source(round(new_x), round(new_y), c);
                    }

                    out(a, b, c) = out(a, b, c) + weight(new_x, new_y) * pixel_value;

                }

            }
        }
    }
}


void applyHomographyMax(const Image &source, const Image &weight, const Image &outweight, Image &out, Matrix &H, bool bilinear) {
    // // --------- HANDOUT  PS07 ------------------------------

    Matrix inv_H = H.inverse();

    cout << "TODO: WEIRD IMAGE BOUNDS THING " << endl;

    cout << "weight image width: " << weight.width() << " height: " << weight.height() << endl;

    cout << "out image width: " << out.width() << " height: " << out.height() << endl;


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

                if (new_x >= 0 && new_y >= 0 && new_x < source.width() - 1 && new_y < source.height() - 1){
                    float pixel_value;
                    if (bilinear){
                        pixel_value = interpolateLin(source, new_x, new_y, c, false);
                        //interpolateLin(const Image &im, float x, float y, int z, bool clamp)
                    }
                    else {
                        pixel_value = source(round(new_x), round(new_y), c);
                    }

                    if (weight(new_x,new_y) > outweight(new_x, new_y)){
                        out(a, b, c) = pixel_value;
                    }

                }

            }
        }
    }
}


Image stitchLinearBlending(const Image &im1, const Image &im2, const Image &we1, const Image &we2, Matrix H) {
    // // --------- HANDOUT  PS07 ------------------------------
    // stitch using image weights.
    // note there is no weight normalization.
    
    BoundingBox B_1 = computeTransformedBBox(im1.width(), im1.height(), H);
    BoundingBox B_2 = computeTransformedBBox(im2.width(), im2.height(), Matrix::Identity(3,3));
    BoundingBox B = bboxUnion(B_1,B_2);
    Matrix T = makeTranslation(B);

    Matrix TH = T*H;

    std::cout << "HELO HUNAM" << endl;
    Image out(B.x2 - B.x1, B.y2 - B.y1, im1.channels());    
    std::cout << "HELO HUNAM 1" << endl;

    //applyhomographyBlend

    applyhomographyBlend(im2, we2, out, T, true);
    std::cout << "HELO HUNAM 2" << endl;

    applyhomographyBlend(im1, we1, out, TH, true);
    std::cout << "HELO HUNAM 2" << endl;

    return out;
}


/*****************************************************************************
 * blending functions Pset 08
 *****************************************************************************/


// low freq and high freq (2-scale decomposition)
vector<Image> scaledecomp(const Image &im, float sigma) {
    vector <Image> ims;
    ims.push_back(gaussianBlur_separable(im, sigma));
    ims.push_back(im - ims[0]);
    return ims;
}

// stitch using different blending models
// blend can be 0 (none), 1 (linear) or 2 (2-layer)
Image stitchBlending(Image &im1, Image &im2, Matrix H, int blend) {
    // // --------- HANDOUT  PS07 ------------------------------
    
    if (blend == 0) {        
        BoundingBox B_1 = computeTransformedBBox(im1.width(), im1.height(), H);
        BoundingBox B_2 = computeTransformedBBox(im2.width(), im2.height(), Matrix::Identity(3,3));
        BoundingBox B = bboxUnion(B_1,B_2);
        Matrix T = makeTranslation(B);

        Image out(B.x2 - B.x1, B.y2 - B.y1, im1.channels());    
        applyHomographyFast(im2, T, out, true);
        applyHomographyFast(im1, T*H, out, true);
        return out;
    }
    else if (blend == 1){

        Image we1 = blendingweight(im1.width(), im1.height());
        Image we2 = blendingweight(im2.width(), im2.height());
        
        BoundingBox B_1 = computeTransformedBBox(im1.width(), im1.height(), H);
        BoundingBox B_2 = computeTransformedBBox(im2.width(), im2.height(), Matrix::Identity(3,3));
        BoundingBox B = bboxUnion(B_1,B_2);
        Matrix T = makeTranslation(B);

        Matrix TH = T*H;

        Image stiched_linear_blending = stitchLinearBlending(im1, im2, we1, we2, H);


        Image out_weight(B.x2 - B.x1, B.y2 - B.y1, 1);
        Image im1_ones(im1.width(), im1.height(), im1.channels());
        im1_ones = 1 - im1_ones;
        Image im2_ones(im2.width(), im2.height(), im2.channels());
        im2_ones = 1 - im2_ones;
        
        applyhomographyBlend(im2_ones, we2, out_weight, T, true);
        applyhomographyBlend(im1_ones, we1, out_weight, TH, true);

        Image out_weight_factor = out_weight;

        for (int i = 0; i < out_weight_factor.width(); i++){
            for (int j = 0; j < out_weight_factor.height(); j++){
                for (int c = 0; c < out_weight_factor.channels(); c++){
                    if (out_weight_factor(i,j,c) == 0) {
                        out_weight_factor(i,j,c) = 1;
                    }
                    else {
                        out_weight_factor(i,j,c) = 1.0/out_weight_factor(i,j,c);
                    }
                }
            }
        }

        Image out = stiched_linear_blending;
        
        for (int i = 0; i < stiched_linear_blending.width(); i++){
            for (int j = 0; j < stiched_linear_blending.height(); j++){
                for (int c = 0; c < stiched_linear_blending.channels(); c++){
                    out(i,j,c) = stiched_linear_blending(i,j,c)*out_weight_factor(i,j);
                }
            }
        }

        return out;
    }
    else if (blend == 2){
        //First, decompose each source image into 
        //low frequencies and high frequencies using 
        //a Gaussian blur. Use a spatial sigma of 2 pixels.
        
        Image we1 = blendingweight(im1.width(), im1.height());
        Image we2 = blendingweight(im2.width(), im2.height());

        float spacial_sigma = 2;

        Image im1_lowfreq = gaussianBlur_separable(im1, spacial_sigma);
        Image im2_lowfreq = gaussianBlur_separable(im2, spacial_sigma);

        Image im1_highfreq = im1 - im1_lowfreq;
        Image im2_highfreq = im2 - im2_lowfreq;

        // For the low frequencies, use the same transition as above.

        Image low_freq_blended_image = stitchBlending(im1_lowfreq, im2_lowfreq, H, 1);

        //For the high frequencies, use an abrupt transition that
        //only keeps the high frequency of the image with the highest weight.

        BoundingBox B_1 = computeTransformedBBox(im1.width(), im1.height(), H);
        BoundingBox B_2 = computeTransformedBBox(im2.width(), im2.height(), Matrix::Identity(3,3));
        BoundingBox B = bboxUnion(B_1,B_2);
        Matrix T = makeTranslation(B);
        Matrix TH = T*H;

        Image high_freq_blended_image(B.x2 - B.x1, B.y2 - B.y1, im1.channels());    
        applyHomographyFast(im2_highfreq, T, high_freq_blended_image, true);
        //applyHomographyMax(im1_highfreq, we1, we2, high_freq_blended_image, TH, true);

        //applyhomographyMax(im1_highfreq, T*H, high_freq_blended_image, true);

        //Compute the final image by adding the resulting low and high fre-quencies.
        Image out = low_freq_blended_image + high_freq_blended_image;
        return out;

    }
    else {
        std::cout << "Error, invalid value of blend";
        return im1;
    }
}

// auto stitch
Image autostitch(Image &im1, Image &im2, int blend, float blurDescriptor, float radiusDescriptor) {
    // // --------- HANDOUT  PS07 ------------------------------
    vector<Point> corners_1 = HarrisCorners(im1);
    vector<Point> corners_2 = HarrisCorners(im2);
    vector<Feature> features_1 = computeFeatures(im1, corners_1, blurDescriptor, radiusDescriptor);
    vector<Feature> features_2 = computeFeatures(im2, corners_2, blurDescriptor, radiusDescriptor);
    vector <FeatureCorrespondence> listOfFeatureCorrespondences = findCorrespondences(features_1, features_2);
    Matrix H = RANSAC(listOfFeatureCorrespondences);

    Image out = stitchBlending(im1, im2, H, blend);

    return out;
}

/************************************************************************
 * Tricks: mini planets.
 ************************************************************************/

vector<float> rect_coords_new_to_old_2(Image old_image, float input_x, float input_y, int newImSize){
    vector <float> radius_and_angle;

    float new_radius = pow((input_x - newImSize/2)*(input_x - newImSize/2) + (input_y - newImSize/2)*(input_y - newImSize/2), .5);
    float new_angle = atan2((input_y - newImSize/2),fabs(input_x - newImSize/2));
    float signed_new_angle = (input_x - newImSize/2) > 0 ? new_angle : -1*new_angle;

    float old_image_radius = pow(old_image.width()*old_image.width() + old_image.height()*old_image.height(),.5);
    float new_img_radius = pow(2*(newImSize/2)*(newImSize/2), .5);

    float radius = new_radius*old_image_radius/new_img_radius;
    float angle = signed_new_angle;

    float old_x = radius*cos(angle);
    float old_y = radius*sin(angle);

    radius_and_angle.push_back(old_x);
    radius_and_angle.push_back(old_y);
    
    return radius_and_angle;
}


Image pano2planet(const Image &pano, int newImSize, bool clamp) {
    // // --------- HANDOUT  PS07 ------------------------------
    
    /*
    Implement Image pano2planet(const Image &pano,
    int newImSize, bool clamp=true).  Map the bottom
    of your input panorama to the center of the the new image and the top
    to a radius corresponding to the distance between the center and the
    right edge (in the square output).
    */

    /*
    Make a new image
    of square size(newImSize),and for each pixel(x, y)in the new image, 
    compute the polar coordinates (angle, radius) assuming that the center
    is the floating point center as in blendingweights.
    */

    cout << "in pano2planet beginning" << endl;

    Image new_img(newImSize, newImSize, pano.channels());

    /*
    Map the bottom of your input panorama to the center of the the new 
    image and the top to a radius corresponding to the distance between
    the center and the right edge (in the square output).

    The left and right sides of the input panorama should be mapped to
    an angle of 0, along the right horizontal axis in the new image with
    increasing (counter-clockwise) angle in the output corresponding to
    sweeping from left to right of the input panorama. Assume standard
    polar coordinate conventions (angle is 0 along right horizontal axis
    and π is along the top vertical axis). Use interpolateLin to copy 2
    pixels from panorama to planet image. Hint: see C++’s atan2.
    */

    for (int i = 0; i < new_img.width(); i++){
        for (int j = 0; j < new_img.height(); j++){
            for (int c = 0; c < new_img.channels(); c++){
                vector<float> new_x_new_y = rect_coords_new_to_old_2(pano, i, j, newImSize);
                float new_x = new_x_new_y[0];
                float new_y = new_x_new_y[1];
                if (new_x >= 0 && new_y >= 0 && new_x < pano.width() - 1 && new_y < pano.height() - 1){
                    float pixel_value = interpolateLin(pano, new_x, new_y, c, false);
                    new_img(i,j,c) = pixel_value;
                }
            }
        }
    }

    cout << "in pano2planet end" << endl;


    return new_img;
}


/************************************************************************
 * 6.865: Stitch N images into a panorama
 ************************************************************************/

// Pset08-865. Compute sequence of N-1 homographies going from Im_i to Im_{i+1}
// Implement me!
vector<Matrix> sequenceHs(vector<Image> ims, float blurDescriptor, float radiusDescriptor) {
    // // --------- HANDOUT  PS07 ------------------------------
    vector<Matrix> Hs;

    //Computes a sequence of N-1 homographies for N images. H[i] should take ims[i] to ims[i+1]. 

    for (int i = 0; i < ims.size() - 1; i++){

        Image im1 = ims[i];
        Image im2 = ims[i];            

        vector<Point> corners_1 = HarrisCorners(im1);
        vector<Point> corners_2 = HarrisCorners(im2);

        vector<Feature> features_1 = computeFeatures(im1, corners_1, blurDescriptor, radiusDescriptor);
        vector<Feature> features_2 = computeFeatures(im2, corners_2, blurDescriptor, radiusDescriptor);

        vector <FeatureCorrespondence> listOfFeatureCorrespondences = findCorrespondences(features_1, features_2);

        Matrix H = RANSAC(listOfFeatureCorrespondences);
        Hs.push_back(H);
    }
    
    return Hs;
}

// stack homographies:
//   transform a list of (N-1) homographies that go from I_i to I_i+1
//   to a list of N homographies going from I_i to I_refIndex.
vector <Matrix> stackHomographies(vector <Matrix> Hs, int refIndex) {
    // // --------- HANDOUT  PS07 ------------------------------
    vector<Matrix> stackedHomographies;
    //Idea: we can create the stacked homographies (going from image i to I_i to I_refIndex) by stacking to that point
    //So, if it's index 3, need to multiply H_0, H_1, H_2, and H_3
    //How can we do this efficiently? Start at the refIndex (whose H is the identity?)
    //Consider: for the sequenceHs, the first H is the identity (kinda)
    //so...we can push them back, but they'll be out of order. Maybe make an n element array?
    //Not going to do the nice way for right now, just the brute forcy way
    
    //Phase 1: index up from 0 to refIndex
    cout << "In stackHomographies. refIndex is: " << refIndex << endl;

    for (int i = 0; i < refIndex; i++){
        Matrix H = Matrix::Identity(3, 3);
        for (int j = i; j < refIndex; j++){
            cout << "In multiplication step of stack homographies" << endl;
            H = H*Hs[j];
        }
        stackedHomographies.push_back(H);
    }

    cout << "Pushed back the indentity matrix at H" << endl;
    //Need to push back the identity as the value at refIndex
    stackedHomographies.push_back(Matrix::Identity(3, 3));

    //Then, all the matrices afterwards
    for (int i = refIndex; i < Hs.size(); i++){
        Matrix H = Matrix::Identity(3, 3);
        cout << "i is: " << i << " and refIndex is: " << refIndex << endl;
        cout << "j should be counting down from i to refIndex" << endl;
        for (int j = i; j >= refIndex; j--){
            cout << "j is: " << j << endl;
            H = Hs[j]*H;
        }
        stackedHomographies.push_back(H);
    }

    cout << "length of stackedHomographies is: " << stackedHomographies.size() << " should be: " << Hs.size() + 1 << endl;
    return stackedHomographies;
}


// Pset08-865: compute bbox around N images given one main reference.
BoundingBox bboxN(const vector<Matrix> &Hs, const vector<Image> &ims) {
    // // --------- HANDOUT  PS07 ------------------------------
    cout << "In bboxN" << endl;
    
    BoundingBox incremental_bbox = computeTransformedBBox(ims[0].width(), ims[0].height(), Hs[0]);
    
    cout << "After incremental_bbox" << endl;

    for (int i = 1; i < ims.size(); i++){
        BoundingBox currentbbox = computeTransformedBBox(ims[i].width(), ims[i].height(), Hs[i]);
        incremental_bbox = bboxUnion(incremental_bbox, currentbbox);
    }

    cout << "After incremental_bbox for loop" << endl;

    
    return incremental_bbox;
}

// Pset08-865.
// Implement me!
Image autostitchN(vector<Image> ims, int refIndex, float blurDescriptor, float radiusDescriptor) {
    // // --------- HANDOUT  PS07 ------------------------------
    /*
    Write autostitch NN, which computes the sequence of homographies using sequenceHs, then propagates
    those homographies using stackHomographies, then computes the overall bounding box, then the 
    translation of the bounding box to (0,0), and finally applies the homographies to all images to get the panorama.
    Use linear blending.
    */

    cout << "Beginning of autostitchNN" << endl;

    cout << "Length of ims is: " << ims.size() << endl;

    vector<Matrix> sequencedHs = sequenceHs(ims, blurDescriptor, radiusDescriptor);

    cout << "After sequencedHs" << endl;

    cout << "Length of sequencedHs is: " << sequencedHs.size() << " should be:  " << ims.size() - 1 << endl;

    vector<Matrix> stackedHomographies = stackHomographies(sequencedHs, refIndex);
    
    cout << "After stackedHomographies" << endl;

    BoundingBox B = bboxN(stackedHomographies, ims);

    cout << "After B" << endl;

    Matrix T = makeTranslation(B);

    cout << "After T" << endl;

    cout << "After first block of autostitchNN" << endl;

    Image ref_image = ims[refIndex];
    Image ref_image_weight = blendingweight(ref_image.width(), ref_image.height());
    Image ref_image_ones(ref_image.width(), ref_image.height(), ref_image.channels());
    ref_image_ones = 1 - ref_image_ones;

    //for out_weight
    Image out_weight_base(B.x2 - B.x1, B.y2 - B.y1, 1);
    applyhomographyBlend(ref_image_ones, ref_image_weight, out_weight_base, T, true);

    Image out(B.x2 - B.x1, B.y2 - B.y1, ref_image.channels());    

    cout << "Made it past autostitchNN setup" << endl;

    for (int i = 0; i < ims.size(); i++){
        Image current_image = ims[i];
        Image current_image_weight = blendingweight(current_image.width(), current_image.height());
        Matrix H = stackedHomographies[i];
        Matrix TH = T*H;
        Image stiched_linear_blending = stitchLinearBlending(current_image, ref_image, current_image_weight, ref_image_weight, TH);

        Image current_image_ones(current_image.width(), current_image.height(), current_image.channels());
        current_image_ones = 1 - current_image_ones;

        //Start with the blended weight of the base image        
        Image current_out_weight = out_weight_base;
        applyhomographyBlend(current_image_ones, current_image_weight, current_out_weight, TH, true);

        Image out_weight_factor = current_out_weight;
        for (int i = 0; i < current_out_weight.width(); i++){
            for (int j = 0; j < current_out_weight.height(); j++){
                for (int c = 0; c < current_out_weight.channels(); c++){
                    if (current_out_weight(i,j,c) == 0) {
                        out_weight_factor(i,j,c) = 1;
                    }
                    else {
                        out_weight_factor(i,j,c) = 1.0/out_weight_factor(i,j,c);
                    }
                }
            }
        }

        out = stiched_linear_blending;
        
        for (int i = 0; i < stiched_linear_blending.width(); i++){
            for (int j = 0; j < stiched_linear_blending.height(); j++){
                for (int c = 0; c < stiched_linear_blending.channels(); c++){
                    out(i,j,c) = stiched_linear_blending(i,j,c)*out_weight_factor(i,j);
                }
            }
        }
        
    }

    return out;

}


/******************************************************************************
 * Helper functions
 *****************************************************************************/

// copy a single-channeled image to several channels
Image copychannels(const Image &im, int nChannels) {
    // image must have one channel
    assert(im.channels() == 1);
    Image oim(im.width(), im.height(), nChannels);

    for (int i = 0; i < im.width(); i++) {
        for (int j = 0; j < im.height(); j++) {
            for (int c = 0; c < nChannels; c++) {
                oim(i, j, c) = im(i, j);
            }
        }
    }
    return oim;
}

