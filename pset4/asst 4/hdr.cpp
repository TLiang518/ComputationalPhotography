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
    for (int i = 0; i < imSeq.size(); i++){
        float mini = epsilonMini;
        if (i == imSeq.size() - 1){
            mini = FLT_MIN;
        }
        float maxi = epsilonMaxi;
        if (i == 0){
            maxi = FLT_MAX;
        }
        weightSeq.push_back(computeWeight(imSeq[i], mini, maxi));
    }

    vector<float> weight_factors_individual;
    vector<float> weight_factors_cumulative;
    float cumulative_factor = 1.0;
    weight_factors_cumulative.push_back(cumulative_factor);
    for (int a = 0; a < weightSeq.size() - 1; a++){    
        float factor = computeFactor(imSeq[a], weightSeq[a], imSeq[a+1], weightSeq[a+1]);
        weight_factors_individual.push_back(factor);
        
        cumulative_factor *= factor;
        weight_factors_cumulative.push_back(cumulative_factor);
    }
    
    Image darkest = imSeq[0];
    Image output(imSeq[0].width(), imSeq[0].height(), imSeq[0].channels());
    for (int a = 0; a < darkest.width(); a++){
        for (int b = 0; b < darkest.height(); b++){
            for (int c = 0; c < darkest.channels(); c++){
                float pixel_value_sum = 0;
                float valid_counter = 0;
                for (int i = 0; i < imSeq.size(); i++){
                    float pixel_weight = weightSeq[i](a,b,c);
                    float pixel_value = imSeq[i](a,b,c);
                    if (pixel_weight > 0){
                        pixel_value_sum += pixel_value/weight_factors_cumulative[i];
                        valid_counter++;
                    }
                }
                if (valid_counter == 0){
                    output(a,b,c) = darkest(a,b,c);
                }
                else {
                    output(a,b,c) = pixel_value_sum/valid_counter;
                }
            }
        }
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

    vector<Image> lumi_chromi_vect = lumiChromi(im);
    Image log10_lumi = log10Image(lumi_chromi_vect[0]);
    float sigma = max(im.width(), im.height())/50.0;
    Image blurred_log_lumi = log10_lumi;
    if (useBila){
        //perform bilateral blurring on the image
        //Image bilateral(const Image &im, float sigmaRange, float sigmaDomain, float truncateDomain, bool clamp)
        blurred_log_lumi = bilateral(log10_lumi, sigmaRange, sigma, 3.0);
    }
    else {
        //perform gaussian blurring on the image
        //Image gaussianBlur_separable(const Image &im, float sigma, float truncate, bool clamp){
        std::cout << "In tone map before separable gaussian blur" << std::endl;
        blurred_log_lumi = gaussianBlur_separable(log10_lumi, sigma, 3.0);
    }
    std::cout << "In tone map after separable gaussian blur" << std::endl;
    
    float k = log10(targetBase)/(blurred_log_lumi.max() - blurred_log_lumi.min());
    Image log_targetBaseImage = k*blurred_log_lumi;
    float log_targetBaseMax = log_targetBaseImage.max();
    float targetBaseMax = pow(10, log_targetBaseMax);
    log_targetBaseImage = log_targetBaseImage - log_targetBaseImage.max();

    Image inLogDetail = log10_lumi - blurred_log_lumi;
    Image log_ampedDetails = detailAmp*inLogDetail;

    Image outlog = log_targetBaseImage + log_ampedDetails;
    Image lumi_output_not_normalized = exp10Image(outlog);
    Image output = lumiChromi2rgb(vector<Image>{lumi_output_not_normalized, lumi_chromi_vect[1]});
    return output;
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
    float image_non_zero_min = image_minnonzero(im);
    Image logImage = im;
    for (int i = 0; i < im.number_of_elements(); i++){
        if (im(i) == 0.0){
            logImage(i) = log10(image_non_zero_min);
        }
        else {
            logImage(i) = log10(im(i));                    
        }        
    }
    return logImage;
}

// Image --> 10^Image
Image exp10Image(const Image &im) {
    // --------- HANDOUT  PS04 ------------------------------
    // take an image in log10 domain and transform it back to linear domain.
    // see pow(a, b)
    Image expImage = im;
    for (int i = 0; i < im.number_of_elements(); i++){
        expImage(i) = pow(10.0, im(i));
    }
    return expImage;
}

// min non-zero pixel value of image
float image_minnonzero(const Image &im) {
    // --------- HANDOUT  PS04 ------------------------------
    // return the smallest value in the image that is non-zeros (across all
    // channels too)
    float min = FLT_MAX;
    for (int i = 0; i < im.number_of_elements(); i++){
        if (im(i) < min) {
            min = im(i);
        }
    }
    return min;
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
