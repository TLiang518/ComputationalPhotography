/* --------------------------------------------------------------------------
 * File:    align.cpp
 * Created: 2015-10-01
 * --------------------------------------------------------------------------
 * 
 * 
 * 
 * ------------------------------------------------------------------------*/



#include "align.h"

using namespace std;

Image denoiseSeq(const vector<Image> &imSeq){
    // // --------- HANDOUT  PS03 ------------------------------
    // Basic denoising by computing the average of a sequence of images
    Image output(imSeq.at(0).width(), imSeq.at(0).height(), imSeq.at(0).channels());
    
    for (auto & im : imSeq) {
        output = output + im/imSeq.size();
    }

    return output;
    
}

Image variance(const vector<Image> &imSeq){
    Image orig_expected_value = denoiseSeq(imSeq);

    std::cout << "In variance_squared method" << std::endl;

    vector<Image> inner_expression;
    for (auto & im : imSeq){
        Image difference = orig_expected_value - im;
        inner_expression.push_back(difference * difference);
    }

    Image output(imSeq.at(0).width(), imSeq.at(0).height(), imSeq.at(0).channels());
    for (auto & im: inner_expression){
        output = output + im;
    }

    output = output/(inner_expression.size() - 1);

    for (int i = 0; i < output.number_of_elements(); i++){
        if (output(i) == 0) {
            std::cout << "ZERO" << std::endl;
            output(i) = 0.000000001;
        }
    }

    return output;
}

Image logSNR(const vector<Image> &imSeq, float scale){
    // // --------- HANDOUT  PS03 ------------------------------
    // returns an image visualizing the per-pixel and
    // per-channel log of the signal-to-noise ratio scaled by scale.
    

    std::cout << "In logSNR" << std::endl;

    Image sigma_squared = variance(imSeq);
    sigma_squared.write("./Output/snr_sigma_squared.png");
    
    std::cout << "computed sigma_squared" << std::endl;

    vector<Image> imSeq_squared;
    for (auto & im : imSeq){
        imSeq_squared.push_back(im*im);
    }
    Image expected_im_squared = denoiseSeq(imSeq_squared); 
    expected_im_squared.write("./Output/snr_expected_im_squared.png");

    std::cout << "computed expected_im_squared" << std::endl;

    Image SNR(imSeq.at(0).width(), imSeq.at(0).height(), imSeq.at(0).channels()); 
    for (int a = 0; a < sigma_squared.width(); a++){
        for (int b = 0; b < sigma_squared.height(); b++){
            for (int c = 0; c < sigma_squared.channels(); c++){
                SNR(a,b,c) = 10*log10(expected_im_squared(a,b,c)/sigma_squared(a,b,c));
            }
        }
    }

    SNR = SNR*scale;

    std::cout << "computed SNR" << std::endl;

    return SNR;
}


vector<int> align(const Image &im1, const Image &im2, int maxOffset){
    // // --------- HANDOUT  PS03 ------------------------------
    // returns the (x,y) offset that best aligns im2 to match im1.
    /*
    Write a function vector<int> align(const Image &im1, const Image &im2, int maxOffset=20) 
    in align.cpp that returns the [x, y] offset that best aligns im2 to match im1. 
    Ignore the difference for all the pixels less than or equal to MaxOffset away from the edges.
    
    Use a brute force approach that tries every possible integer translation and evaluates the quality
    of a match using the squared error norm (the sum of the squared pixel differences).
    
    The Image roll(const Image &im, int xRoll, int yRoll) func- tion in align.cpp might come in handy. 
    It circularly shifts an image, causing borders to wrap around. However, since you will be ignoring boundary pixels,
    wrapping the pixel values should not be a problem. Make sure to test your procedure before moving on.
    */

    int x = 0;
    int y = 0;
    float best_squared_error_norm = FLT_MAX;
    for (int x_offset = -maxOffset; x_offset <= maxOffset; x_offset++){
        for (int y_offset = -maxOffset; y_offset <= maxOffset; y_offset++){
            Image rolled_im2 = roll(im2, x_offset, y_offset);
            float squared_error_norm = 0;
            for (int a = maxOffset; a < im1.width() - maxOffset; a++){
                for (int b = maxOffset; b < im1.height() - maxOffset; b++){
                    for (int c = 0; c < im1.channels(); c++){
                        squared_error_norm += pow(rolled_im2(a,b,c) - im1(a,b,c),2);
                    }
                }
            }
            if (squared_error_norm < best_squared_error_norm) {
                best_squared_error_norm = squared_error_norm;
                x = x_offset;
                y = y_offset; 
            }
        }


    }
    return vector<int> {x,y};
}

Image alignAndDenoise(const vector<Image> &imSeq, int maxOffset){
    // // --------- HANDOUT  PS03 ------------------------------
    // Registers all images to the first one in a sequence and outputs
    // a denoised image even when the input sequence is not perfectly aligned.
    Image im1 = imSeq.at(0);

    vector<Image> shifted_seq;
    for (auto & im : imSeq){
        std::vector<int> best_shift = align(im1, im);
        Image shifted_image = roll(im, best_shift.at(0), best_shift.at(1));
        shifted_seq.push_back(shifted_image);
    }
    
    Image denoised_image = denoiseSeq(shifted_seq);
    return denoised_image;
}

Image split(const Image &sergeyImg){
    // --------- HANDOUT  PS03 ------------------------------
    // 6.865 only:
    // split a Sergey images to turn it into one 3-channel image.
    Image reconstructedSergey(sergeyImg.width(), floor(sergeyImg.height()/3.0), 3);
    float new_height = floor(sergeyImg.height()/3.0);
    for (int a = 0; a < sergeyImg.width(); a++){
        for (int b = 0; b < new_height*3; b++){
            if (b < new_height){
                reconstructedSergey(a,b,0) = sergeyImg(a,b);
            }
            else if (b >= new_height && b < 2*new_height){
                reconstructedSergey(a,b - new_height, 1) = sergeyImg(a,b);
            }
            else if (b >= new_height) {
                reconstructedSergey(a,b - 2*new_height, 2) = sergeyImg(a,b);
            }
            else {
                std::cout << "New height is: " << new_height << std::endl;
                std::cout << "Failed to fall into right case. Indicies: " << a << " " << b << std::endl;
            }
        }
    }
    return reconstructedSergey;
}

Image sergeyRGB(const Image &sergeyImg, int maxOffset){
    // // --------- HANDOUT  PS03 ------------------------------
    // 6.865 only:
    // aligns the green and blue channels of your rgb channel of a sergey
    // image to the red channel. This should return an aligned RGB image
    Image splitSergey = split(sergeyImg);
    Image sergeyRed(splitSergey.width(), splitSergey.height(), 1);    
    Image sergeyGreen(splitSergey.width(), splitSergey.height(), 1);
    Image sergeyBlue(splitSergey.width(), splitSergey.height(), 1);
    
    for (int a = 0; a < splitSergey.width(); a++){
        for (int b = 0; b < splitSergey.height(); b++){
            sergeyRed(a,b) = splitSergey(a,b,0);
            sergeyGreen(a,b) = splitSergey(a,b,1);
            sergeyBlue(a,b) = splitSergey(a,b,2);
        }
    }
    vector<int> green_align = align(sergeyRed, sergeyGreen, maxOffset);
    vector<int> blue_align = align(sergeyRed, sergeyBlue, maxOffset);

    Image rolled_sergeyGreen = roll(sergeyGreen, green_align[0], green_align[1]);
    Image rolled_sergeyBlue = roll(sergeyBlue, blue_align[0], blue_align[1]);

    Image sergeyAligned = splitSergey;
    for (int a = maxOffset; a < splitSergey.width() - maxOffset; a++){
        for (int b = maxOffset; b < splitSergey.height() - maxOffset; b++){
            sergeyAligned(a,b,0) = sergeyRed(a,b);
            sergeyAligned(a,b,1) = rolled_sergeyGreen(a, b);
            sergeyAligned(a,b,2) = rolled_sergeyBlue(a, b);
        }
    }
    return sergeyAligned;

}


/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// This circularly shifts an image by xRoll in the x direction and
// yRoll in the y direction. xRoll and yRoll can be negative or postive
Image roll(const Image &im, int xRoll, int yRoll){
    
    int xNew, yNew;
    Image imRoll(im.width(), im.height(), im.channels());
    
    // for each pixel in the original image find where its corresponding
    // location is in the rolled image
    for (int x=0; x<im.width(); x++){
        for (int y=0; y<im.height(); y++){
            
            // use modulo to figure out where the new location is in the
            // rolled image. Then take care of when this returns a negative number
            xNew = (x + xRoll) % im.width();
            yNew = (y + yRoll) % im.height();
            xNew = (xNew<0)*(imRoll.width() + xNew) + (xNew>=0)*xNew;
            yNew = (yNew<0)*(imRoll.height() + yNew) + (yNew>=0)*yNew;
            
            // assign the rgb values for each pixel in the original image to
            // the location in the new image
            for (int z=0; z<im.channels(); z++){
                imRoll(xNew, yNew, z) = im(x,y,z);
            }
        }
    }
    
    // return the rolled image
    return imRoll;
}
