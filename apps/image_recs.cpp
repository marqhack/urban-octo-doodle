/**
 *  Copyright 2015 Mike Reed
 */

#include "image.h"
#include "GCanvas.h"
#include "GBitmap.h"
#include "GColor.h"
#include "GPoint.h"
#include "GRect.h"
#include "GShader.h"
#include <string>

#define RAMP_W      1
#define RAMP_H     28

static void draw_solid_ramp(GCanvas* canvas) {
    const float c = 1.0 / 512;
    const float d = 1.0 / 256;

    const struct {
        GColor  fC0, fDC;
    } rec[] = {
        { GColor::MakeARGB(1,   c,   c,   c), GColor::MakeARGB(0,  d,  d,  d) },   // grey
        { GColor::MakeARGB(1, 1-c,   0,   0), GColor::MakeARGB(0, -d,  0,  0) },   // red
        { GColor::MakeARGB(1,   0,   c,   c), GColor::MakeARGB(0,  0,  d,  d) },   // cyan
        { GColor::MakeARGB(1,   0, 1-c,   0), GColor::MakeARGB(0,  0, -d,  0) },   // green
        { GColor::MakeARGB(1,   c,   0,   c), GColor::MakeARGB(0,  d,  0,  d) },   // magenta
        { GColor::MakeARGB(1,   0,   0, 1-c), GColor::MakeARGB(0,  0,  0, -d) },   // blue
        { GColor::MakeARGB(1,   c,   c,   0), GColor::MakeARGB(0,  d,  d,  0) },   // yellow
    };

    
    for (int y = 0; y < GARRAY_COUNT(rec); ++y) {
        GColor color = rec[y].fC0;
        GColor delta = rec[y].fDC;
        for (int x = 0; x < 256; x++) {
            const GRect rect = GRect::MakeXYWH(x * RAMP_W, y * RAMP_H, RAMP_W, RAMP_H);
            canvas->fillRect(rect, color);
            color.fA += delta.fA;
            color.fR += delta.fR;
            color.fG += delta.fG;
            color.fB += delta.fB;
        }
    }
}

static void offset(GRect* r, float dx, float dy) {
    r->fLeft += dx;
    r->fRight += dx;
    r->fTop += dy;
    r->fBottom += dy;
}

static void draw_blend_ramp(GCanvas* canvas, const GColor& bg) {
    canvas->clear(bg);

    GRect rect = GRect::MakeXYWH(-25, -25, 70, 70);

    int delta = 8;
    for (int i = 0; i < 200; i += delta) {
        float r = i / 200.0;
        float g = fabs(cos(i/40.0));
        float b = fabs(sin(i/50.0));
        GColor color = GColor::MakeARGB(0.3, r, g, b);
        canvas->fillRect(rect, color);
        offset(&rect, delta, delta);
    }
}

static void draw_blend_white(GCanvas* canvas) {
    draw_blend_ramp(canvas, GColor::MakeARGB(1, 1, 1, 1));
}

static void draw_blend_black(GCanvas* canvas) {
    draw_blend_ramp(canvas, GColor::MakeARGB(1, 0, 0, 0));
}

///////////////////////////////////////////////////////////////////////////////////////////////////

static void draw_spocks_quad(GCanvas* canvas) {
    const int N = 300;

    GBitmap tex;
    tex.readFromFile("apps/spock.png");

    for (int y = 0; y < 2; ++y) {
        for (int x = 0; x < 2; ++x) {
            canvas->fillBitmapRect(tex, GRect::MakeXYWH(x * N - N/2, y * N - N/2, N, N));
        }
    }
}

static void draw_spocks_zoom(GCanvas* canvas) {
    const int N = 300;
    
    GBitmap tex;
    tex.readFromFile("apps/spock.png");

    for (int i = 0; i < 9; ++i) {
        GRect r = GRect::MakeLTRB(i * 10, i * 10, N - i * 10, N - i * 10);
        canvas->fillBitmapRect(tex, r);
    }
}

// After scaling by this, the caller need just cast to (int)
static const float gScaleUnitToByte = 255.99999f;

static GPixel pin_and_premul_to_pixel(GColor c) {
    c = c.pinToUnit();
    
    float a = c.fA * gScaleUnitToByte;
    int ia = (int)a;
    int ir = (int)(a * c.fR);
    int ig = (int)(a * c.fG);
    int ib = (int)(a * c.fB);
    return GPixel_PackARGB(ia, ir, ig, ib);
}

static void make_circle(const GBitmap& bitmap, const GColor& color) {
    const GPixel px = pin_and_premul_to_pixel(color);
    
    const float cx = (float)bitmap.width() / 2;
    const float cy = (float)bitmap.height() / 2;
    const float radius = cx - 1;
    const float radius2 = radius * radius;
    
    GPixel* dst = bitmap.pixels();
    for (int y = 0; y < bitmap.height(); ++y) {
        const float dy = y - cy;
        for (int x = 0; x < bitmap.width(); ++x) {
            const float dx = x - cx;
            const float dist2 = dx*dx + dy*dy;
            if (dist2 <= radius2) {
                dst[x] = px;
            } else {
                dst[x] = 0; // transparent
            }
        }
        dst = (GPixel*)((char*)dst + bitmap.rowBytes());
    }
}

class AutoBitmap : public GBitmap {
public:
    AutoBitmap(int width, int height) {
        // just to exercise the ability to have a rowbytes > width
        const int slop = (height >> 1) * sizeof(GPixel);

        fWidth = width;
        fHeight = height;
        fRowBytes = fWidth * sizeof(GPixel) + slop;
        fPixels = (GPixel*)malloc(fRowBytes * fHeight);
    }

    ~AutoBitmap() {
        free(fPixels);
    }
};

static void draw_bm_circles(GCanvas* canvas) {
    const int N = 300;

    AutoBitmap src(N, N);

    const struct {
        GRect   fRect;
        GColor  fColor;
    } recs[] = {
        { GRect::MakeXYWH(  0,   0,   N,   N), GColor::MakeARGB(1, 1, 1, 1) },

        { GRect::MakeXYWH(  0,   0, N/2, N/2), GColor::MakeARGB(0.8f, 0, 0, 1) },
        { GRect::MakeXYWH(N/2,   0, N/2, N/2), GColor::MakeARGB(0.6f, 0, 1, 0) },
        { GRect::MakeXYWH(  0, N/2, N/2, N/2), GColor::MakeARGB(0.4f, 1, 0, 0) },
        { GRect::MakeXYWH(N/2, N/2, N/2, N/2), GColor::MakeARGB(0.2f, 0, 0, 0) },

        { GRect::MakeXYWH(  0, N/3,   N, N/3), GColor::MakeARGB(0.5f, 1, 1, 0) },
        { GRect::MakeXYWH(N/3,   0, N/3,   N), GColor::MakeARGB(0.5f, 0, 1, 1) },
        { GRect::MakeXYWH(N/3, N/3, N/3, N/3), GColor::MakeARGB(0.5f, 1, 0, 1) },
    };

    for (int i = 0; i < GARRAY_COUNT(recs); ++i) {
        make_circle(src, recs[i].fColor);
        canvas->fillBitmapRect(src, recs[i].fRect);
    }
}

static void draw_circle_big(GCanvas* canvas) {
    const int N = 300;
    const float alpha = 0.4f;
    const GColor colors[] = {
        GColor::MakeARGB(alpha, 1, 0, 0),
        GColor::MakeARGB(alpha, 0, 1, 0),
        GColor::MakeARGB(alpha, 0, 0, 1),
    };

    int x = 0;
    int n = N;
    for (int i = 0; n > 4; ++i) {
        AutoBitmap src(n, n);
        make_circle(src, colors[i % GARRAY_COUNT(colors)]);
        canvas->fillBitmapRect(src, GRect::MakeXYWH(x, 0, N, N));
        x += N / 12;
        n >>= 1;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

static void draw_tri(GCanvas* canvas) {
    GPoint pts[] = {
        { 10, 10 },
        { 200, 50 },
        { 100, 200 },
    };
    canvas->fillConvexPolygon(pts, GARRAY_COUNT(pts), GColor::MakeARGB(1, 0, 1, 0));
}

static void draw_tri_clipped(GCanvas* canvas) {
    GPoint pts[] = {
        { -10, -10 },
        { 300, 50 },
        { 100, 300 },
    };
    canvas->fillConvexPolygon(pts, GARRAY_COUNT(pts), GColor::MakeARGB(1, 1, 1, 0));
}

static void make_regular_poly(GPoint pts[], int count, float cx, float cy, float radius) {
    float angle = 0;
    const float deltaAngle = M_PI * 2 / count;
    
    for (int i = 0; i < count; ++i) {
        pts[i].set(cx + cos(angle) * radius, cy + sin(angle) * radius);
        angle += deltaAngle;
    }
}

static void dr_poly(GCanvas* canvas, float dx, float dy) {
    GPoint storage[12];
    for (int count = 12; count >= 3; --count) {
        make_regular_poly(storage, count, 256, 256, count * 10 + 120);
        for (int i = 0; i < count; ++i) {
            storage[i].fX += dx;
            storage[i].fY += dy;
        }
        GColor c = GColor::MakeARGB(0.8f,
                                    fabs(sin(count*7)),
                                    fabs(sin(count*11)),
                                    fabs(sin(count*17)));
        canvas->fillConvexPolygon(storage, count, c);
    }
}

static void draw_poly(GCanvas* canvas) {
    dr_poly(canvas, 0, 0);
}

static void draw_poly_center(GCanvas* canvas) {
    dr_poly(canvas, -128, -128);
}

static GPoint scale(GPoint vec, float size) {
    float scale = size / sqrt(vec.fX * vec.fX + vec.fY * vec.fY);
    return GPoint::Make(vec.fX * scale, vec.fY * scale);
}

static void draw_line(GCanvas* canvas, GPoint a, GPoint b, float width, const GColor& color) {
    GPoint norm = scale(GPoint::Make(b.fY - a.fY, a.fX - b.fX), width/2);
    
    GPoint pts[4];
    pts[0] = GPoint::Make(a.fX + norm.fX, a.fY + norm.fY);
    pts[1] = GPoint::Make(b.fX + norm.fX, b.fY + norm.fY);
    pts[2] = GPoint::Make(b.fX - norm.fX, b.fY - norm.fY);
    pts[3] = GPoint::Make(a.fX - norm.fX, a.fY - norm.fY);

    canvas->fillConvexPolygon(pts, 4, color);
}

static void draw_poly_rotate(GCanvas* canvas) {
    const GPoint start = GPoint::Make(20, 20);
    const float scale = 200;

    const int N = 10;
    GColor color = GColor::MakeARGB(1, 1, 0, 0);
    const float deltaR = -1.0 / N;
    const float deltaB = 1.0 / N;
    
    const float width = 10;

    for (float angle = 0; angle <= M_PI/2; angle += M_PI/2/N) {
        GPoint end = GPoint::Make(start.fX + cos(angle) * scale,
                                  start.fY + sin(angle) * scale);
        draw_line(canvas, start, end, width, color);

        color.fR += deltaR;
        color.fB += deltaB;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

#define DRAW_CONCAT_W   800
#define DRAW_CONCAT_H   600

/*
 *  A(1 - t)^2 + 2Bt(1 - t) + Ct^2
 *
 *  A - 2tA + At^2 + 2Bt - 2Bt^2 + Ct^2
 *
 *  t^0 : A
 *  t^1 : -2A + 2B
 *  t^2 : A - 2B + C
 */
static void parametric2coeff(float p0, float p1, float p2, float coeff[3]) {
    coeff[0] = p0;
    coeff[1] = 2 * (p1 - p0);
    coeff[2] = p0 - 2 * p1 + p2;
}

/**
 *  f(t) = At^2 + Bt + C
 *       = (At + B)t + C
 */
static float eval_quad(const float coeff[3], float t) {
    return (coeff[2] * t + coeff[1]) * t + coeff[0];
}
static GPoint eval_quad(const GPoint quad[3], float t) {
    float coeffX[3], coeffY[3];

    parametric2coeff(quad[0].fX, quad[1].fX, quad[2].fX, coeffX);
    parametric2coeff(quad[0].fY, quad[1].fY, quad[2].fY, coeffY);

    return GPoint::Make(eval_quad(coeffX, t), eval_quad(coeffY, t));
}

static float interp_color(float v0, float v1, float t) {
    return v0 + (v1 - v0) * t;
}

static GColor interp_color(const GColor colors[2], float t) {
    float radians = t * 4 * M_PI;
    float value = cos(radians);
    t = (value + 1) * 0.5f;
    return GColor::MakeARGB(interp_color(colors[0].fA, colors[1].fA, t),
                            interp_color(colors[0].fR, colors[1].fR, t),
                            interp_color(colors[0].fG, colors[1].fG, t),
                            interp_color(colors[0].fB, colors[1].fB, t));
}

class Drawable {
public:
    virtual void draw(GCanvas*, float t) const = 0;
};

template <typename T> void draw_concat(GCanvas* canvas, T xform, const Drawable& drawable) {
    const GPoint quad[] = {
        { 100, 500 }, { 400, -200 }, { 700, 500 }
    };
    for (float t = 0; t <= 1; t += (1.0f/32)) {
        const GPoint pt = eval_quad(quad, t);
        canvas->save();
        canvas->translate(pt.fX, pt.fY);
        xform(canvas, t);
        drawable.draw(canvas, t);
        canvas->restore();
    }
}

class ColorsDrawable : public Drawable {
    const GColor* fColors;
    GPoint fPts[20];
    int fCount;

public:
    ColorsDrawable(const GColor colors[], int count) : fColors(colors) {
        GASSERT(count <= GARRAY_COUNT(fPts));
        make_regular_poly(fPts, count, 0, 0, DRAW_CONCAT_W/5);
        fCount = count;
    }
    void draw(GCanvas* canvas, float t) const override {
        canvas->fillConvexPolygon(fPts, fCount, interp_color(fColors, t));
    }
};

class BitmapDrawable : public Drawable {
    GBitmap fBitmap;
public:
    BitmapDrawable() {
        fBitmap.readFromFile("apps/spock.png");
    }
    void draw(GCanvas* canvas, float) const override {
        const GRect rect = GRect::MakeLTRB(-DRAW_CONCAT_W/10, -DRAW_CONCAT_H/10,
                                           DRAW_CONCAT_W/10,  DRAW_CONCAT_H/10);
        canvas->fillBitmapRect(fBitmap, rect);
    }
};

static void draw_concat_scale(GCanvas* canvas) {
    const GColor colors[] = {
        GColor::MakeARGB(1, 0, 0, 1), GColor::MakeARGB(1, 0, 1, 0)
    };
    draw_concat(canvas,
                [](GCanvas* canvas, float t) {
                    canvas->scale(0.5f + t, 0.5f + t*(1 - t));
                }, ColorsDrawable(colors, 7));
}

static void draw_concat_scale_bitmap(GCanvas* canvas) {
    draw_concat(canvas,
                [](GCanvas* canvas, float t) {
                    canvas->scale(0.5f + t, 0.5f + t*(1 - t));
                }, BitmapDrawable());
}

static void draw_concat_rotate(GCanvas* canvas) {
    const GColor colors[] = {
        GColor::MakeARGB(1,0, 1, 0), GColor::MakeARGB(1, 1, 0, 0)
    };
    draw_concat(canvas,
                [](GCanvas* canvas, float t) {
                    canvas->rotate(t * 4 * M_PI);
                }, ColorsDrawable(colors, 5));
}

static void draw_concat_rotate_bitmap(GCanvas* canvas) {
    draw_concat(canvas,
                [](GCanvas* canvas, float t) {
                    canvas->rotate(t * 4 * M_PI);
                }, BitmapDrawable());
}

///////////////////////////////////////////////////////////////////////////////////////////////////

static void draw_grad_insets(GCanvas* canvas) {
    const GColor colors[] = {
        GColor::MakeARGB(1, 1, 0, 0), GColor::MakeARGB(1, 0, 0, 1)
    };
    
    GRect r = GRect::MakeXYWH(50, 10, 300, 70);

    for (int inset = -100; inset <= 100; inset += 50) {
        const GPoint pts[] = {{ r.left() + inset, 0 }, { r.right() - inset, 0 }};
        GShader* shader = GShader::FromLinearGradient(pts, colors);
        canvas->shadeRect(r, shader);
        delete shader;
        r.offset(0, r.height() + 2);
    }
}

static void draw_grad_insets_rot(GCanvas* canvas) {
    canvas->translate(400, 0);
    canvas->rotate(M_PI / 2);
    draw_grad_insets(canvas);
}

static void draw_bm_localmatrix(GCanvas* canvas, const GRect& r, const GBitmap& tex,
                                const float localMatrix[6]) {
    GShader* shader = GShader::FromBitmap(tex, localMatrix);
    canvas->shadeRect(r, shader);
    delete shader;
}

static void draw_bm_pad(GCanvas* canvas) {
    GBitmap tex;
    tex.readFromFile("apps/spock.png");

    const float scale = 0.5f;
    const float w = tex.width() * scale;
    const float h = tex.height() * scale;
    const float extra = 100;

    const float mat0[] = {
        scale, 0, 10,
        0, scale, 10
    };
    draw_bm_localmatrix(canvas, GRect::MakeXYWH(10, 10, w, h), tex, mat0);

    canvas->save();
    canvas->translate(w + 20, 0);
    const float mat1[] = {
        scale, 0, 10 + extra/2,
        0, scale, 10
    };
    draw_bm_localmatrix(canvas, GRect::MakeXYWH(10, 10, w + extra, h), tex, mat1);
    canvas->restore();
    
    canvas->save();
    canvas->translate(0, h + 20);
    const float mat2[] = {
        scale, 0, 10,
        0, scale, 10 + extra/2,
    };
    draw_bm_localmatrix(canvas, GRect::MakeXYWH(10, 10, w, h + extra), tex, mat2);
    canvas->restore();
    
    canvas->save();
    canvas->translate(w + 20, h + 20);
    const float mat3[] = {
        scale, 0, 10 + extra/2,
        0, scale, 10 + extra/2,
    };
    draw_bm_localmatrix(canvas, GRect::MakeXYWH(10, 10, w + extra, h + extra), tex, mat3);
    canvas->restore();
}

static void draw_regular_poly(GCanvas* canvas, int count, float cx, float cy, float radius,
                              GShader* shader) {
    GPoint* pts = new GPoint[count];
    make_regular_poly(pts, count, cx, cy, radius);
    canvas->shadeConvexPolygon(pts, count, shader);
    delete[] pts;
    delete shader;
}

static void draw_poly_shader(GCanvas* canvas, GShader* shader) {
    draw_regular_poly(canvas, 7, 100, 100, 90, shader);
}

static void draw_poly_shaders(GCanvas* canvas) {
    GBitmap tex;
    tex.readFromFile("apps/spock.png");

    for (int repeat = 0; repeat < 2; ++repeat) {
        const float mat[] = {
            0.75f, 0, -50,
            0, 0.75f, -20,
        };
        draw_poly_shader(canvas, GShader::FromBitmap(tex, mat));

        canvas->translate(200, 0);

        const GColor colors[] = {
            GColor::MakeARGB(1, 0, 1, 0), GColor::MakeARGB(1, 0, 0, 1)
        };
        const GPoint pts[] = {
            GPoint::Make(0, 200), GPoint::Make(200, 0)
        };
        draw_poly_shader(canvas, GShader::FromLinearGradient(pts, colors));

        canvas->translate(200, 400);
        canvas->scale(-1, -1);
    }
}

static void draw_radial(GCanvas* canvas) {
    const GColor colors[] = {
        GColor::MakeARGB(1, 1, 0, 1), GColor::MakeARGB(1, 0.2, 0.5, 0.2)
    };
    GShader* shader = GShader::FromRadialGradient(GPoint::Make(200, 200), 150, colors);
    draw_regular_poly(canvas, 128, 200, 200, 150, shader);
}

static void draw_radial_quad(GCanvas* canvas) {
    const GColor colors[] = {
        GColor::MakeARGB(1, 0.9, 1, 0.9), GColor::MakeARGB(1, 0.1, 0.5, 0.1)
    };
    GShader* shader = GShader::FromRadialGradient(GPoint::Make(0, 0), 50, colors);
    const GRect r = GRect::MakeLTRB(-50, -50, 50, 50);
    
    canvas->scale(4, 4);
    
    canvas->shadeRect(r, shader);
    canvas->translate(100, 0);
    canvas->shadeRect(r, shader);
    canvas->translate(-100, 100);
    canvas->shadeRect(r, shader);
    canvas->translate(100, 0);
    canvas->shadeRect(r, shader);
    
    delete shader;
}

static void draw_linear_big(GCanvas* canvas) {
    const GColor colors[] = {
        GColor::MakeARGB(1, 1, 0, 0), GColor::MakeARGB(1, 0, 0, 1)
    };
    const GPoint pts[] = {
        GPoint::Make(150, 180), GPoint::Make(250, 220)
    };
    GShader* shader = GShader::FromLinearGradient(pts, colors);
    canvas->shadeRect(GRect::MakeWH(400, 400), shader);
    delete shader;
}

static void draw_radial_big(GCanvas* canvas) {
    const GColor colors[] = {
        GColor::MakeARGB(1, 1, 0, 0), GColor::MakeARGB(1, 0, 0, 1)
    };
    GShader* shader = GShader::FromRadialGradient(GPoint::Make(200, 200), 70, colors);
    canvas->shadeRect(GRect::MakeWH(400, 400), shader);
    delete shader;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

static void show_line(GCanvas* canvas, const GPoint& p0, const GPoint& p1,
                      const GCanvas::Stroke& stroke, const GColor& color) {
    canvas->strokeLine(p0, p1, stroke, color);

    GColor c = GColor::MakeARGB(1, 0, 1, 0);
    canvas->fillRect(GRect::MakeXYWH(p0.x() - 2, p0.y() - 2, 5, 5), c);
    canvas->fillRect(GRect::MakeXYWH(p1.x() - 2, p1.y() - 2, 5, 5), c);
}

static void draw_stroke_lines(GCanvas* canvas) {
    GCanvas::Stroke stroke{ 50, 0, false };
    const GPoint p0 = GPoint::Make(80, 200);
    const GPoint p1 = GPoint::Make(300, 80);

    show_line(canvas, p0, p1, stroke, GColor::MakeARGB(1, 1, 0, 0));
    canvas->translate(50, 100);
    stroke.fAddCap = true;
    show_line(canvas, p1, p0, stroke, GColor::MakeARGB(1, 0, 0, 1));
}

static float function(float x) {
    return sin(x);// / x;
}

static void draw_stroke_polyline(GCanvas* canvas) {
//    return;

    const int N = 100;
    GPoint pts[N];
    float x = -2*M_PI;
    float dx = -2*x / (N - 1);
    for (int i = 0; i < N; ++i) {
        pts[i] = GPoint::Make(x, function(x));
        x += dx;
    }

    for (float angle = 0; angle < M_PI; angle += M_PI / 6) {
        GShader* shader = GShader::FromColor(GColor::MakeARGB(1, angle/4, sin(angle), sin(angle+2)));
        canvas->save();
        canvas->translate(200, 200);
        canvas->scale(30, 30);
        canvas->rotate(angle);
        canvas->strokePolygon(pts, N, false, GCanvas::Stroke{0.25f,0,true}, shader);
        canvas->restore();
    }
}

static void draw_stroke_rects(GCanvas* canvas) {
    const GRect r = GRect::MakeXYWH(50, 60, 100, 80);

    canvas->strokeRect(r, GCanvas::Stroke{30, 2, false}, GColor::MakeARGB(1, 1, 0, 0));
    canvas->translate(200, 0);
    canvas->strokeRect(r, GCanvas::Stroke{60, 2, false}, GColor::MakeARGB(1, 0, 1, 0));
    canvas->translate(-200, 200);
    canvas->strokeRect(r, GCanvas::Stroke{60, 1, false}, GColor::MakeARGB(1, 0, 0, 1));
    canvas->translate(200, 0);
    canvas->strokeRect(r, GCanvas::Stroke{30, 1, false}, GColor::MakeARGB(1, 0, 0, 0));
}

static void draw_strokes(GCanvas* canvas) {
    GShader* shader = GShader::FromColor(GColor::MakeARGB(1, 1, 0, 0));
    GCanvas::Stroke stroke = { 20, 2, false };

    canvas->strokeRect(GRect::MakeXYWH(30, 30, 300, 200), stroke, shader);
}

static void make_regular_star(GPoint pts[], int count, float cx, float cy, float radius, bool cw) {
    const int n = count >> 1;

    float angle = -M_PI / 2;
    float deltaAngle = 2 * M_PI * n / count;
    if (!cw) {
        deltaAngle = -deltaAngle;
    }

    for (int i = 0; i < count; ++i) {
        pts[i].set(cx + cos(angle) * radius, cy + sin(angle) * radius);
        angle += deltaAngle;
    }
}

static void draw_star(GCanvas* canvas, const GColor& color, float miter, bool cw) {
    const int n = 7;
    GPoint pts[n];
    const GCanvas::Stroke stroke = { 20, miter, false };

    make_regular_star(pts, n, 200, 200, 140, cw);
    GShader* shader = GShader::FromColor(color);
    canvas->strokePolygon(pts, n, true, stroke, shader);
    delete shader;
}

static void draw_stroke_star(GCanvas* canvas) {
    const struct {
        GPoint  fOffset;
        GColor  fColor;
        float   fMiterLimit;
        bool    fClockwise;
    } recs[] = {
        { GPoint::Make(  0,   0), GColor::MakeARGB(1, 1, 0, 0), 8, true },
        { GPoint::Make(400,   0), GColor::MakeARGB(1, 0, 1, 0), 8, false },
        { GPoint::Make(  0, 400), GColor::MakeARGB(1, 0, 0, 1), 1, true },
        { GPoint::Make(400, 400), GColor::MakeARGB(1, 0, 0, 0), 1, false },
    };

    canvas->scale(0.5f, 0.5f);
    for (const auto& r : recs) {
        canvas->save();
        canvas->translate(r.fOffset.x(), r.fOffset.y());
        draw_star(canvas, r.fColor, r.fMiterLimit, r.fClockwise);
        canvas->restore();
    }
}

static void shade_stroke_ring(GCanvas* canvas, GShader* shader, float width) {
    canvas->save();
    canvas->scale(1, -1);
    canvas->translate(0, -400);

    const int n = 8;
    GPoint poly[n];
    const GCanvas::Stroke stroke = { width, 3, true };
    
    make_regular_poly(poly, n, 200, 200, 150);
    canvas->strokePolygon(poly, n, true, stroke, shader);

    canvas->restore();
}

static GShader* make_shader0() {
    const GColor colors[] = {
        GColor::MakeARGB(1, 1, 0, 0), GColor::MakeARGB(1, 0, 0, 1)
    };
    const GPoint pts[] = { { 0, 0 }, { 400, 400 } };
    return GShader::FromLinearGradient(pts, colors);
}

static GShader* make_shader1() {
    const GColor colors[] = {
        GColor::MakeARGB(1, 1, 1, 1), GColor::MakeARGB(1, 0, 0, 0)
    };
    return GShader::FromRadialGradient(GPoint::Make(200, 200), 250, colors);
}

static GShader* make_shader2() {
    GBitmap tex;
    tex.readFromFile("apps/spock.png");
    return GShader::FromBitmap(tex, GRect::MakeWH(400, 400));
}

static void draw_stroke_rings(GCanvas* canvas) {
    const struct {
        GShader* (*fMakeShader)();
        float    fWidth;
        GPoint   fOffset;
    } recs[] = {
        { make_shader0, 80, GPoint::Make(0, 0), },
        { make_shader1, 80, GPoint::Make(400, 0), },
        { make_shader2, 80, GPoint::Make(200, 400), },
    };

    canvas->scale(0.5f, 0.5f);
    for (const auto& r : recs) {
        canvas->save();
        canvas->translate(r.fOffset.x(), r.fOffset.y());
        GShader* shader = r.fMakeShader();
        shade_stroke_ring(canvas, shader, r.fWidth);
        delete shader;
        canvas->restore();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

const CS575DrawRec gDrawRecs[] = {
    { draw_solid_ramp,  256 * RAMP_W, 7*RAMP_H, 1, "solid_ramp"    },
    { draw_blend_white, 200, 200,               1, "blend_white"   },
    { draw_blend_black, 200, 200,               1, "blend_black"   },

    { draw_spocks_quad, 300, 300,               2, "spocks_quad"   },
    { draw_spocks_zoom, 300, 300,               2, "spocks_zoom"   },
    { draw_bm_circles,  300, 300,               2, "circles_blend" },
    { draw_circle_big,  400, 300,               2, "circles_fat"   },

    { draw_tri,         256, 256,               3, "tri"           },
    { draw_tri_clipped, 256, 256,               3, "tri_clipped"   },
    { draw_poly,        512, 512,               3, "poly"          },
    { draw_poly_center, 256, 256,               3, "poly_center"   },
    { draw_poly_rotate, 230, 230,               3, "poly_rotate"   },

    { draw_concat_scale,        DRAW_CONCAT_W, DRAW_CONCAT_H, 4, "draw_concat_scale" },
    { draw_concat_scale_bitmap, DRAW_CONCAT_W, DRAW_CONCAT_H, 4, "draw_concat_scale_bitmap" },
    { draw_concat_rotate,       DRAW_CONCAT_W, DRAW_CONCAT_H, 4, "draw_concat_rotate" },
    { draw_concat_rotate_bitmap,DRAW_CONCAT_W, DRAW_CONCAT_H, 4, "draw_concat_rotate_bitmap" },

    { draw_grad_insets, 400, 400,               5, "gradient_insets" },
    { draw_grad_insets_rot, 400, 400,           5, "gradient_insets_rot" },
    { draw_bm_pad, 530, 500,                    5, "bitmap_shader_clamp" },
    { draw_poly_shaders, 400, 400,              5, "poly_shaders" },
    { draw_radial,      400, 400,               5, "radial" },
    { draw_radial_quad, 400, 400,               5, "radial_quad" },
    { draw_linear_big,  400, 400,               5, "linear_big" },
    { draw_radial_big,  400, 400,               5, "radial_big" },

    { draw_stroke_lines,     400, 400,          6, "stroke_line" },
    { draw_stroke_polyline,  400, 400,          6, "stroke_polyline" },
    { draw_stroke_rects,     400, 400,          6, "stroke_rect" },
    { draw_stroke_star,      400, 400,          6, "stroke_star" },
    { draw_stroke_rings,     400, 400,          6, "stroke_rings" },

    { NULL, 0, 0, 0, NULL },
};
