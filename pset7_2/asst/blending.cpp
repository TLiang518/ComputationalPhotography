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


}

Image stitchLinearBlending(const Image &im1, const Image &im2, const Image &we1, const Image &we2, Matrix H) {
    // // --------- HANDOUT  PS07 ------------------------------
    // stitch using image weights.
    // note there is no weight normalization.
    return Image(1,1,1);
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
    return Image(1,1,1);
}

// auto stitch
Image autostitch(Image &im1, Image &im2, int blend, float blurDescriptor, float radiusDescriptor) {
    // // --------- HANDOUT  PS07 ------------------------------
    return Image(1,1,1);
}

/************************************************************************
 * Tricks: mini planets.
 ************************************************************************/

Image pano2planet(const Image &pano, int newImSize, bool clamp) {
    // // --------- HANDOUT  PS07 ------------------------------
    return Image(1,1,1);
}


/************************************************************************
 * 6.865: Stitch N images into a panorama
 ************************************************************************/

// Pset08-865. Compute sequence of N-1 homographies going from Im_i to Im_{i+1}
// Implement me!
vector<Matrix> sequenceHs(vector<Image> ims, float blurDescriptor, float radiusDescriptor) {
    // // --------- HANDOUT  PS07 ------------------------------
    return vector<Matrix>();
}

// stack homographies:
//   transform a list of (N-1) homographies that go from I_i to I_i+1
//   to a list of N homographies going from I_i to I_refIndex.
vector <Matrix> stackHomographies(vector <Matrix> Hs, int refIndex) {
    // // --------- HANDOUT  PS07 ------------------------------
    return vector<Matrix>();
}


// Pset08-865: compute bbox around N images given one main reference.
BoundingBox bboxN(const vector<Matrix> &Hs, const vector<Image> &ims) {
    // // --------- HANDOUT  PS07 ------------------------------
    return BoundingBox(0,0,0,0);
}

// Pset08-865.
// Implement me!
Image autostitchN(vector<Image> ims, int refIndex, float blurDescriptor, float radiusDescriptor) {
    // // --------- HANDOUT  PS07 ------------------------------
    return Image(1,1,1);
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

