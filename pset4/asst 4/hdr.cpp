// hdr.cpp
// Assignment 5


#include "hdr.h"
#include "filtering.h"
#include <math.h>
#include <algorithm>


using namespace std;

/**************************************************************
 //                       HDR MERGING                        //
 *************************************************************/

Image computeWeight(const Image &im, float epsilonMini, float epsilonMaxi){
    // --------- HANDOUT  PS04 ------------------------------
    // Generate a weight image that indicates which pixels are good to use in
    // HDR, i.e. weight=1 when the pixel value is in [epsilonMini, epsilonMaxi].
    // The weight is per pixel, per channel.
    Image output = im;

    for (int a = 0; a < im.width(); a++){
        for (int b = 0; b < im.height(); b++){
            for (int c = 0; c < im.channels(); c++){
                if (im(a,b,c) >= epsilonMini && im(a,b,c) <= epsilonMaxi) {
                    output(a,b,c) = 1.0;
                }
                else {
                    output(a,b,c) = 0.0;
                }
            }
        }
    }
    return output;
}


float computeFactor(const Image &im1, const Image &w1, const Image &im2, const Image &w2){
    // --------- HANDOUT  PS04 ------------------------------
    // Compute the multiplication factor between a pair of images. This
    // gives us the relative exposure between im1 and im2. It is computed as 
    // the median of im2/(im1+eps) for some small eps, taking into account
    // pixels that are valid in both images.
    float eps = 10e-10;
    vector<float> w1_factors;
    vector<float> w2_factors;

    for (int a = 0; a < im1.width(); a++){
        for (int b = 0; b < im1.height(); b++){
            for (int c = 0; c < im1.channels(); c++){
                if (w1(a,b,c) == 1.0 && w2(a,b,c) == 1.0){
                    w1_factors.push_back(im1(a,b,c) + eps);
                    w2_factors.push_back(im2(a,b,c) + eps);
                }
            }
        }
    }
    sort(w1_factors.begin(), w1_factors.end());
    sort(w2_factors.begin(), w2_factors.end());
    int index = floor(w2_factors.size()/2.0);
    return w2_factors[index]/w1_factors[index];
}


Image makeHDR(vector<Image> &imSeq, float epsilonMini, float epsilonMaxi){
    // --------- HANDOUT  PS04 ------------------------------
    // Merge images to make a single hdr image
    // For each image in the sequence, compute the weight map (special cases
    // for the first and last images).
    // Compute the exposure factor for each consecutive pair of image.
    // Write the valid pixel to your hdr output, taking care of rescaling them
    // properly using the factor.
    vector<Image> weightSeq;

    for (auto & im : imSeq){
        weightSeq.push_back(computeWeight(im));
    }

    vector<Image> normalized_weightSeq;
    normalized_weightSeq.push_back(weightSeq[0]);
    for (int a = 1; a < weightSeq.size() - 1; a++){    
        float exposure_factor_difference = computeFactor(imSeq[a], weightSeq[a], imSeq[a+1], weightSeq[a+1]);
        Image normalized_image = imSeq[a]/exposure_factor_difference;
        normalized_weightSeq.push_back(normalized_image);
    }
    normalized_weightSeq.push_back(weightSeq[0]);

    for (int a = 0; a < normalized_weightSeq.size(); a++){
        normalized_weightSeq[a].write("./Output/test_" + to_string(a));
    }

    Image output(imSeq[0].width(), imSeq[0].height(), imSeq[0].channels());
    int counter = 0;
    for (auto & normalized_im : normalized_weightSeq) {
        output = output + normalized_im;
        counter++;
        output.write("./Output/output_" + to_string(counter));
    }
    return output;
}

/**************************************************************
 //                      TONE MAPPING                        //
 *************************************************************/


Image toneMap(const Image &im, float targetBase, float detailAmp, bool useBila, float sigmaRange) {
    // --------- HANDOUT  PS04 ------------------------------
    // tone map an hdr image
    // - Split the image into its luminance-chrominance components.
    // - Work in the log10 domain for the luminance
    // - 
    return im;

}



/*********************************************************************
 *                       Tone mapping helpers                        *
 *********************************************************************/



// image --> log10Image
Image log10Image(const Image &im) {
    // --------- HANDOUT  PS04 ------------------------------
    // Taking a linear image im, transform to log10 scale.
    // To avoid infinity issues, make any 0-valued pixel be equal the the minimum
    // non-zero value. See image_minnonzero(im).
    return im;
    
}

// Image --> 10^Image
Image exp10Image(const Image &im) {
    // --------- HANDOUT  PS04 ------------------------------
    // take an image in log10 domain and transform it back to linear domain.
    // see pow(a, b)
    return im;
    
}

// min non-zero pixel value of image
float image_minnonzero(const Image &im) {
    // --------- HANDOUT  PS04 ------------------------------
    // return the smallest value in the image that is non-zeros (across all
    // channels too)
    return 0.0f;
    
}

/*********************************************************************
 *                       Do not edit below                           *
 *********************************************************************/
Image changeGamma(const Image & im, float old_gamma, float new_gamma) {
    // Image output(im.width(), im.height(), im.channels());
    // Figure out what power to take the values of im, to get the values of output
    // return output;
    float exponent = new_gamma/old_gamma;
    Image output = im;
    for (int i = 0 ; i < im.number_of_elements();i++) {
        output(i) = pow(im(i), exponent);
    }
    return output;
}

