/*
 *  Copyright 2015 Mike Reed
 */

#ifndef GCanvas_DEFINED
#define GCanvas_DEFINED

#include "GTypes.h"

class GBitmap;
class GColor;
class GPoint;
class GRect;
class GShader;

class GCanvas {
public:
    static GCanvas* Create(const GBitmap&);

    virtual ~GCanvas() {}

    /**
     *  Saves a copy of the CTM, allowing subsequent modifications (by calling concat()) to be
     *  undone when restore() is called.
     *
     *  e.g.
     *  // CTM is in state A
     *  canvas->save();
     *  canvas->conact(...);
     *  canvas->conact(...);
     *  // Now the CTM has been modified by the calls to concat()
     *  canvas->restore();
     *  // Now the CTM is again in state A
     */
    virtual void save() = 0;

    /**
     *  Balances calls to save(), returning the CTM to the state it was in when the corresponding
     *  call to save() was made. These calls can be nested.
     *
     *  e.g.
     *  save()
     *      concat()
     *      concat()
     *      save()
     *          concat()
     *          ...
     *      restore()
     *      ...
     *  restore()
     */
    virtual void restore() = 0;

    /**
     *  Modifies the CTM (current transformation matrix) by pre-concatenating it with the specfied
     *  matrix.
     *
     *  CTM' = CTM * matrix
     *
     *  The result is that any drawing that uses the new CTM' will be affected AS-IF it were
     *  first transformed by matrix, and then transformed by the previous CTM.
     */
    virtual void concat(const float matrix[6]) = 0;

    /**
     *  Pretranslates the CTM by the specified tx, ty
     */
    void translate(float tx, float ty);

    /**
     *  Prescales the CTM by the specified sx, sy
     */
    void scale(float sx, float sy);

    /**
     *  Prerotates the CTM by the specified angle in radians.
     */
    void rotate(float radians);

    /**
     *  Fill the entire canvas with the specified color.
     *
     *  This completely overwrites any previous colors, it does not blend.
     */
    virtual void clear(const GColor&) = 0;
    
    /**
     *  Fill the rectangle with the color.
     *
     *  The affected pixels are those whose centers are "contained" inside the rectangle:
     *      e.g. contained == center > min_edge && center <= max_edge
     *
     *  Any area in the rectangle that is outside of the bounds of the canvas is ignored.
     *
     *  If the color's alpha is < 1, blend it using SRCOVER blend mode.
     */
    virtual void fillRect(const GRect&, const GColor&) = 0;

    /**
     *  Scale and translate the bitmap such that it fills the specific rectangle.
     *
     *  Any area in the rectangle that is outside of the bounds of the canvas is ignored.
     *
     *  If a given pixel in the bitmap is not opaque (e.g. GPixel_GetA() < 255) then blend it
     *  using SRCOVER blend mode.
     */
    virtual void fillBitmapRect(const GBitmap& src, const GRect& dst) = 0;

    /**
     *  Fill the convex polygon with the color, following the same "containment" rule as
     *  rectangles.
     *
     *  Any area in the polygon that is outside of the bounds of the canvas is ignored.
     *
     *  If the color's alpha is < 1, blend it using SRCOVER blend mode.
     */
    virtual void fillConvexPolygon(const GPoint[], int count, const GColor&) = 0;

    /**
     *  Fill the specified rect using the shader. The colors returned by the shader are blended
     *  into the canvas using SRC_OVER blend mode.
     */
    virtual void shadeRect(const GRect& rect, GShader* shader) = 0;

    /**
     *  Fill the specified polygon using the shader. The colors returned by the shader are blended
     *  into the canvas using SRC_OVER blend mode.
     */
    virtual void shadeConvexPolygon(const GPoint[], int count, GShader* shader) = 0;

    struct Stroke {
        float   fWidth;         // width of the stroke
        float   fMiterLimit;    // limiting ratio of diagonal length of a miterjoin to strokewidth.
                                // If that ratio exceeds fMiterLimit, then the join should be drawn
                                // with a bevel.
        bool    fAddCap;        // if true, add a square cap at the ends
    };

    void strokeLine(const GPoint& p0, const GPoint& p1, const Stroke&, const GColor&);
    void strokeLine(const GPoint& p0, const GPoint& p1, const Stroke&, GShader*);
    void strokeRect(const GRect&, const Stroke&, const GColor&);
    void strokeRect(const GRect&, const Stroke&, GShader*);

    /**
     *  Stroke the specified polygon using the Stroke settings. If isClosed is true, then the
     *  drawn stroke should connect the first and last points of the polygon, else it should not,
     *  and those end-caps should reflect the Stroke.fAddCap setting.
     */
    virtual void strokePolygon(const GPoint[], int n, bool isClosed, const Stroke&, GShader*) = 0;
};

#endif
