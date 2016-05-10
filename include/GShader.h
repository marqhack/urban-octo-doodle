/*
 *  Copyright 2015 Mike Reed
 */

#ifndef GShader_DEFINED
#define GShader_DEFINED

#include "GPixel.h"

class GBitmap;
class GColor;
class GPoint;
class GRect;

/**
 *  GShaders create colors to fill whatever geometry is being drawn to a GCanvas.
 */
class GShader {
public:
    virtual ~GShader() {}

    /**
     *  Called before each use, this tells the shader the CTM for the current drawing.
     *  This returns true if the shader can handle the CTM, and therefore it is valid to call
     *  shadeRow(). If it cannot handle the CTM, this will return false, and shadeRow()
     *  should not be called.
     */
    virtual bool setContext(const float ctm[6]) = 0;

    /**
     *  Given a row of pixels in device space [x, y] ... [x + count - 1, y], return the
     *  corresponding src pixels in row[0...count - 1]. The caller must ensure that row[]
     *  can hold at least [count] entries.
     */
    virtual void shadeRow(int x, int y, int count, GPixel row[]) = 0;

    /** Return a subclass of GShader that draws the specified color. */
    static GShader* FromColor(const GColor&);

    /**
     *  Return a subclass of GShader that draws the specified bitmap and local-matrix.
     *  Returns null if the either parameter is not valid.
     */
    static GShader* FromBitmap(const GBitmap&, const float localMatrix[6]);
    static GShader* FromBitmap(const GBitmap&, const GRect&);

    /**
     *  Return a subclass of GShader that draws the specified linear gradient. Returns NULL if
     *  the parameters are not valid.
     */
    static GShader* FromLinearGradient(const GPoint pts[2], const GColor colors[2]);
    
    /**
     *  Return a subclass of GShader that draws the specified radial gradient. Returns NULL if
     *  the parameters are not valid.
     */
    static GShader* FromRadialGradient(const GPoint& center, float radius, const GColor colors[2]);
};

#endif
