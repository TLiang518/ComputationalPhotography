#ifndef A10_H_PHUDVTKB
#define A10_H_PHUDVTKB

#include "Image.h"
#include "basicImageManipulation.h"
#include "filtering.h"
#include "matrix.h"
#include <iostream>
#include <cmath>

// Write your declarations here, or extend the Makefile if you add source
// files
void brush(Image & im, int x, int y, vector<float> color, Image &texture);
void singleScalePaint(Image & im, Image & out, Image & importance, Image & texture, int N = 10000, int size = 50, float noise = 0.3);
Image sharpnessMap(Image &im, float sigma = 1.0);
Image anisotropic_gaussian(Image &im, float sigma = 1.0);

Image painterly(Image &im, Image &texture, int N = 10000, int size = 50, float noise = 0.3);

//Harris Corner Response code from pset7
//Original Image computeTensor(const Image &im, float sigmaG=1, float factorSigma=4);
Image computeTensor(const Image &im, float sigmaG=3, float factorSigma=5);
Image computeAngles(Image & im);
void singleScaleOrientedPaint(Image & im, Image & out, Image & importance, Image & texture, int N = 10000, int size = 50, float noise = 0.3, int nAngles = 36, bool useRegularStroke = true);
vector<Image> rotateBrushes(Image &texture, int n = 36);
vector<Image> rotateBrushesPerpendicular(Image &texture, int n = 36);
Image orientedPaint(Image &im, Image &texture, int N = 10000, int size = 50, float noise = 0.3, bool useRegularStroke = true);

#endif /* end of include guard: A10_H_PHUDVTKB */

