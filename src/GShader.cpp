/*
 *  Copyright 2015 Mike Reed
 */

#include "GShader.h"
#include "GBitmap.h"
#include "GColor.h"
#include "GPixel.h"
#include "GRect.h"

namespace {

static GPixel pin_and_premul_to_pixel(GColor c) {
    const float scaleUnitToByte = 255.9999f;
    
    // normalize it first, so all components are "in range"
    c = c.pinToUnit();
    
    float a = c.fA * scaleUnitToByte;
    int ia = (int)a;
    int ir = (int)(a * c.fR);
    int ig = (int)(a * c.fG);
    int ib = (int)(a * c.fB);
    return GPixel_PackARGB(ia, ir, ig, ib);
}

class PixelShader : public GShader {
public:
    PixelShader(GPixel src) : fSrc(src) {}
    
    bool setContext(const float ctm[6]) override {
        // Since we're just a single color, we can ignore the ctm parameter.
        return true;
    }
    
    void shadeRow(int, int, int count, GPixel row[]) override {
        for (int i = 0; i < count; ++i) {
            row[i] = fSrc;
        }
    }
    
private:
    GPixel  fSrc;
};

} // namespace

GShader* GShader::FromColor(const GColor& color) {
    return new PixelShader(pin_and_premul_to_pixel(color));
}

GShader* GShader::FromBitmap(const GBitmap& bm, const GRect& dst) {
    const float mat[] = {
        dst.width() / bm.width(),   0,                              dst.left(),
        0,                          dst.height() / bm.height(),     dst.top(),
    };
    return FromBitmap(bm, mat);
}
