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
    return Vec2f(-a.y, a.x);
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
    lPQ = length(subtracted_vec);

    e1 = scalarMult(subtracted_vec, 1/lPQ);
    //http://gamedev.stackexchange.com/questions/70075/how-can-i-find-the-perpendicular-to-a-2d-vector
    Vec2f e2_intermediate = perpendicular(e1);
    e2 = scalarMult(perpendicular(subtracted_vec), 1/lPQ);

    std::cout << "P " << "(" << P.x << " , " << P.y << ")" << " Q " <<  "(" << Q.x << " , " << Q.y << ")"  << " lPQ " << lPQ << endl;
    std::cout << "e1 " << "(" << e1.x << " , " << e1.y << ")" << " e2 " <<  "(" << e2.x << " , " << e2.y << ")"  << endl;
    std::cout << "Check that e1 and e2 are perpendicular: dot product: " << dot(e1,e2) << endl;
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
    float u = dot(subtract(X,P), subtract(Q, P))/pow(lPQ,2);
    float v = dot(subtract(X,P), perpendicular(subtract(Q, P)))/lPQ;
    std::cout << "In XtoUV, u is: " << u << " and v is: " << v << endl;
    return Vec2f(u,v);
}


Vec2f Segment::UVtoX(Vec2f uv) const {
    // --------- HANDOUT  PS03 ------------------------------
    // compute the (x, y) position of a point given by the (u,v)
    // location relative to this segment.
    // * Be careful with the different normalization for u and v
    float u = uv.x;
    float v = uv.y;
    Vec2f intermediate_X = scalarMult(e1, lPQ*u);
    Vec2f X = add(add(intermediate_X, scalarMult(e2, v)), P);

    cout << "In XtoUV, X is: (" << X.x << "," << X.y << ")" << endl;

    return X;
}


float Segment::distance(Vec2f X) const {
    // // --------- HANDOUT  PS03 ------------------------------
    // // Implement distance from a point X(x,y) to the segment. Remember the 3
    // // cases from class.

    //Three cases -- one, the length of the perpendicular line segment will be the optimal distance
    // U is the vector perpendicular to P
    // Second and third cases -- the endpoints of the line segments are idea

    //Following general example of http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
   
    Vec2f subtracted_vec = subtract(Q, P);
    float length_squared = pow(length(subtracted_vec),2);
    if (length_squared == 0){
        return length(subtract(Q, X));
    }
    //project point p onto line
    float project_p_onto_line = dot(subtract(X,P), subtract(Q,P))/length_squared;

    if (project_p_onto_line < 0.0) {
        return length(subtract(X, P));
    }  
    else if (project_p_onto_line > 1.0) {
        return length(subtract(X, Q));
    } 

    Vec2f projection =  add(P, scalarMult(subtract(Q,P), project_p_onto_line));  // Q + project_p_onto_line * (Q - P)
    return length(subtract(projection, X));
}


Image warpBy1(const Image &im, const Segment &segBefore, const Segment &segAfter){
    // --------- HANDOUT  PS03 ------------------------------
    // Warp an entire image according to a pair of segments.
     // Warp an entire image according to a pair of segments.
    Image output(im.width(), im.height(), im.channels());
    for (int a=0; a < im.width(); a++){
        for (int b=0; b < im.height(); b++){
            Vec2f uv = segAfter.XtoUV(Vec2f(a,b));
            Vec2f X = segBefore.UVtoX(uv);
            for (int c = 0; c < im.channels(); c++){
                output(a,b,c) = interpolateLin(im, X.x, X.y, c, true);
            }
        }
    } 
    return output;
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
