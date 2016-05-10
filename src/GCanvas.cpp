/**
 *  Copyright 2015 Mike Reed
 */

#include "GCanvas.h"
#include "GPoint.h"
#include "GRect.h"
#include "GShader.h"

void GCanvas::translate(float tx, float ty) {
    const float mat[6] = {
        1, 0, tx,
        0, 1, ty,
    };
    this->concat(mat);
}

void GCanvas::scale(float sx, float sy) {
    const float mat[6] = {
        sx, 0, 0,
        0, sy, 0,
    };
    this->concat(mat);
}

void GCanvas::rotate(float radians) {
    const float c = cos(radians);
    const float s = sin(radians);
    const float mat[6] = {
        c, -s, 0,
        s, c, 0,
    };
    this->concat(mat);
}

// All of the methods below eventually call the virtual strokePolygon()

void GCanvas::strokeLine(const GPoint& p0, const GPoint& p1, const Stroke& s, const GColor& color) {
    GShader* shader = GShader::FromColor(color);
    if (shader) {
        this->strokeLine(p0, p1, s, shader);
        delete shader;
    }
}

void GCanvas::strokeLine(const GPoint& p0, const GPoint& p1, const Stroke& s, GShader* shader) {
    const GPoint pts[] = { p0, p1 };
    this->strokePolygon(pts, 2, false, s, shader);
}

void GCanvas::strokeRect(const GRect& rect, const Stroke& stroke, const GColor& color) {
    GShader* shader = GShader::FromColor(color);
    if (shader) {
        this->strokeRect(rect, stroke, shader);
        delete shader;
    }
}

static inline GPoint* to_quad(const GRect& r, GPoint pts[4]) {
    pts[0] = GPoint::Make(r.left(), r.top());
    pts[1] = GPoint::Make(r.right(), r.top());
    pts[2] = GPoint::Make(r.right(), r.bottom());
    pts[3] = GPoint::Make(r.left(), r.bottom());
    return pts;
}

void GCanvas::strokeRect(const GRect& rect, const Stroke& stroke, GShader* shader) {
    GPoint pts[4];
    this->strokePolygon(to_quad(rect, pts), 4, true, stroke, shader);
}
