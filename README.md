# 6.865: Computational Photography (Fall 2015)

The following repository contains my code from 6.865, Computational
Photography, taught by Fredo Durand at MIT in Fall 2015. Below is a quick overview of the 9 problemsets assigned and the 10th problemset (a final project).

## Pset 10: Final Project - Painterly Rendering

For this problemset, we were allowed to choose a topic of our own based off of psets from previous semesters or techniques described in computational photograhy papers. The aim of my project was to take an input photograph and output a version of the photograph styled as if it were painted. Some additions to the project were to include and incorporate edge detection into the algorithm (such that strokes around edges are smaller and denser) and follow edges and orientations within an image by calculating the eigenvalue of the structure tensor to provide a more convincing painterly style. Read my project writeup here: https://github.com/jdesa/ComputationalPhotography/blob/master/pset10/asst/6.865FinalProject.pdf and see the pset10 folder for the code. The "pset<number>.pdf" files in each pset folder contain the assignment information for that pset.

## Psets 8 and 9: Halide

These two problemsets dipped into Halide, a language for efficeint image
processing; Problemset 8 consisted of an introduction to Halide including basic syntax
and scheduling, and Problemset 9 consisted of an implementation of an image processing pipeline of
concepts discussed earlier in the course (such as the Harris corners detector).

## Pset 7 and 7-2: Automatic Panorama Stitching

These combined psets built on the homography code and concepts developed in pset
6 to work on the Harris coner detector, RANSAC, and blending/stitching of images to create a system that can do automatic panorama stitching of multiple images.

## Pset 6: Homographies

This pset served as an introduction to homography, leading up to an
implementation of image stitching that would later lead to the automatic
ponorama stitching assignment in pset 7.

## Pset 5: Warping

This pset involving implementing the warping algorithm described in the Feature
Based Image Metamorphasis Paper (https://github.com/jdesa/ComputationalPhotography/blob/master/pset5/FeatureBasedImageMetamorphosis.pdf).

## Pset 4: HDR and Tonemapping

Merge images to HDR & use tonemapping to display images on low dynamic range
screens.

## Pset 3: Denoising and Demosaicing

Implemented denoising based on averaging of images, variance and signal-to-noise
computation, image alignment using least squares, color channel demosaicing
(green channel, red and blue channel, edge-based green, red and blue channel
based on difference to green).

## Pset 2: Convolution and the Bilateral Filter

Blurring, Image filtering using the Sobel Kernel, Gaussian filtering (1D, 2D,
Separable 2D), Sharpening, Denoising, YUV bilateral gaussian filtering.

## Pset 1: Basic Image Processing

Brightness/Contrast Change, Luminance/Chrominance Decoupling, YUV, White
Balancing of Images, and the Spanish Castle Illusion.