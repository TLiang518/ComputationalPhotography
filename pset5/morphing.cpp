/* -----------------------------------------------------------------
 * File:    morphing.cpp
 * Created: 2015-09-25
 * -----------------------------------------------------------------
 * 
 * 
 * 
 * ---------------------------------------------------------------*/




#include <cassert>
#include "morphing.h"

using namespace std;

Vec2f add(const Vec2f & a, const Vec2f & b) {
    // --------- HANDOUT  PS03 ------------------------------
    // Return the vector sum of a an b
    return Vec2f(a.x + b.x, a.y + b.y);
    
}


Vec2f subtract(const Vec2f & a, const Vec2f & b) {
    // --------- HANDOUT  PS03 ------------------------------
    // Return a-b
    return Vec2f(a.x - b.x, a.y - b.y);
}


Vec2f scalarMult(const Vec2f & a, float f) {
    // --------- HANDOUT  PS03 ------------------------------
    // Return a*f
    return Vec2f(f*a.x,f*a.y);
}


float dot(const Vec2f & a, const Vec2f & b) {
    // --------- HANDOUT  PS03 ------------------------------
    // Return the dot product of a and b.
    return a.x*b.x + a.y*b.y;
}

float length(const Vec2f & a) {
    // --------- HANDOUT  PS03 ------------------------------
    // Return the length of a.
    return sqrt(pow(a.x,2) + pow(a.y,2)); 
}

Vec2f perpendicular(const Vec2f & a){
    return Vec2f(a.y, -a.x);
}


// The Segment constructor takes in 2 points P(x1,y1) and Q(x2,y2) correspoding to
// the ends of a segment and initialize the local reference frame e1,e2.
Segment::Segment(Vec2f P_, Vec2f Q_) : P(P_), Q(Q_) {
    // // --------- HANDOUT  PS03 ------------------------------
    // // The initializer list above ": P(P_), Q(Q_)" already copies P_
    // // and Q_, so you don't have to do it in the body of the constructor.
    // You should:
    // * Initialize the local frame e1,e2 (see header file)
    Vec2f subtracted_vec = subtract(Q, P);
    e1 = Vec2f(subtracted_vec.x/length(subtracted_vec), subtracted_vec.y/length(subtracted_vec));
    //http://gamedev.stackexchange.com/questions/70075/how-can-i-find-the-perpendicular-to-a-2d-vector
    Vec2f e2_intermediate = perpendicular(e1);
    e2 = Vec2f(e2_intermediate.x/length(e2_intermediate), e2_intermediate.y/length(e2_intermediate));
    lPQ = length(subtracted_vec);
}


Vec2f Segment::XtoUV(Vec2f X) const {
    // --------- HANDOUT  PS03 ------------------------------
    // Compute the u,v coordinates of a point X with
    // respect to the local frame of the segment.
    // e2 ^
    //    |
    // v  +  * X
    //    | /    
    //    |/
    //    *--+------>-----*
    //    P  u     e1     Q
    //                    u=1
    //
    // * Be careful with the different normalization for u and v
    
    /*
    Now that our Segment class is usable let’s implement methods to
    convert from the global (x,y) coordinates of a point to the local (u,v)
    coordinates in the reference frame as in Beier and Neely. Be careful:
    although v exactly corresponds to the second coordinate in the local frame,
    u is actually rescaled so that the u coordinate of Q is 1 
    (Equations (1) and (2) in the paper).
    */
    float u = dot(subtract(X,P), subtract(Q, P))/pow(length(subtract(Q, P)),2);
    float v = dot(subtract(X,P), perpendicular(subtract(Q, P)))/length(subtract(Q, P));
    return Vec2f(u,v); 
}


Vec2f Segment::UVtoX(Vec2f uv) const {
    // --------- HANDOUT  PS03 ------------------------------
    // compute the (x, y) position of a point given by the (u,v)
    // location relative to this segment.
    // * Be careful with the different normalization for u and v
    float u = uv.x;
    float v = uv.y;

    float X_x =  u * pow(lPQ,2)/(Q.x - P.x) + P.x;
    float X_y = v * pow(lPQ,2)/(Q.y - P.y) + P.y;

    return Vec2f(X_x, X_y);
}


float Segment::distance(Vec2f X) const {
    // // --------- HANDOUT  PS03 ------------------------------
    // // Implement distance from a point X(x,y) to the segment. Remember the 3
    // // cases from class.

    //Three cases -- one, the length of the perpendicular line segment will be the optimal distance
    // U is the vector perpendicular to P
    // Second and third cases -- the endpoints of the line segments are idea

}


Image warpBy1(const Image &im, const Segment &segBefore, const Segment &segAfter){
    // --------- HANDOUT  PS03 ------------------------------
    // Warp an entire image according to a pair of segments.
    return im;
}


float Segment::weight(Vec2f X, float a, float b, float p) const {
    // --------- HANDOUT  PS03 ------------------------------
    // compute the weight of a segment to a point X(x,y) given the weight
    // parameters a,b, and p (see paper for details).
    return 1.0f; // changeme
}


Image warp(const Image &im, const vector<Segment> &src_segs,
        const vector<Segment> &dst_segs, float a, float b, float p)
{
    // --------- HANDOUT  PS03 ------------------------------
    // Warp an image according to a vector of before and after segments using
    // segment weighting
    return im;
}


vector<Image> morph(const Image &im_before, const Image &im_after,
        const vector<Segment> &segs_before, const vector<Segment> &segs_after, 
        int N, float a, float b, float p)
{
    // --------- HANDOUT  PS03 ------------------------------
    // return a vector of N+2 images: the two inputs plus N images that morphs
    // between im_before and im_after for the corresponding segments. im_before should be the first image, im_after the last.
    return vector<Image>();

}
