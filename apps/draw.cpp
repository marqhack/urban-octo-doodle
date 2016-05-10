/**
 *  Copyright 2015 Mike Reed
 *
 *  COMP 575 -- Fall 2015
 */

#include "GWindow.h"
#include "GBitmap.h"
#include "GCanvas.h"
#include "GColor.h"
#include "GRandom.h"
#include "GRect.h"
#include "GShader.h"

#include <cstdlib>
#include <vector>

static const float CORNER_SIZE = 9;

template <typename T> int find_index(const std::vector<T*>& list, T* target) {
    for (int i = 0; i < list.size(); ++i) {
        if (list[i] == target) {
            return i;
        }
    }
    return -1;
}

static GRandom gRand;

static GColor rand_color() {
    return GColor::MakeARGB(0.5f, gRand.nextF(), gRand.nextF(), gRand.nextF());
}

static GRect make_from_pts(const GPoint& p0, const GPoint& p1) {
    return GRect::MakeLTRB(std::min(p0.fX, p1.fX), std::min(p0.fY, p1.fY),
                           std::max(p0.fX, p1.fX), std::max(p0.fY, p1.fY));
}

static bool contains(const GRect& rect, float x, float y) {
    return rect.left() < x && x < rect.right() && rect.top() < y && y < rect.bottom();
}

static GRect offset(const GRect& rect, float dx, float dy) {
    return GRect::MakeLTRB(rect.left() + dx, rect.top() + dy,
                           rect.right() + dx, rect.bottom() + dy);
}

static bool hit_test(float x0, float y0, float x1, float y1) {
    const float dx = fabs(x1 - x0);
    const float dy = fabs(y1 - y0);
    return std::max(dx, dy) <= CORNER_SIZE;
}

static bool in_resize_corner(const GRect& r, float x, float y, GPoint* anchor) {
    if (hit_test(r.left(), r.top(), x, y)) {
        anchor->set(r.right(), r.bottom());
        return true;
    } else if (hit_test(r.right(), r.top(), x, y)) {
        anchor->set(r.left(), r.bottom());
        return true;
    } else if (hit_test(r.right(), r.bottom(), x, y)) {
        anchor->set(r.left(), r.top());
        return true;
    } else if (hit_test(r.left(), r.bottom(), x, y)) {
        anchor->set(r.right(), r.top());
        return true;
    }
    return false;
}

static void draw_corner(GCanvas* canvas, const GColor& c, float x, float y, float dx, float dy) {
    canvas->fillRect(make_from_pts(GPoint::Make(x, y - 1), GPoint::Make(x + dx, y + 1)), c);
    canvas->fillRect(make_from_pts(GPoint::Make(x - 1, y), GPoint::Make(x + 1, y + dy)), c);
}

class Shape {
public:
    Shape() : fAngle(0) {}
    
    virtual ~Shape() {}
    virtual GRect getRect() = 0;
    virtual void setRect(const GRect&) {}
    virtual GColor getColor() = 0;
    virtual void setColor(const GColor&) {}
    virtual GClick* findClickHandler(GPoint loc) { return nullptr; }
    virtual bool handleKey(int) { return false; }
    virtual void drawHilite(GCanvas* canvas);

    void concat(GCanvas* canvas) {
        const GRect r = this->getRect();
        const float cx = (r.left() + r.right()) / 2;
        const float cy = (r.top() + r.bottom()) / 2;

        canvas->translate(cx, cy);
        canvas->rotate(fAngle * M_PI / 180);
        canvas->translate(-cx, -cy);
    }

    void draw(GCanvas* canvas) {
        canvas->save();
        this->concat(canvas);
        this->onDraw(canvas);
        canvas->restore();
    }
    
    int fAngle;
    
protected:
    virtual void onDraw(GCanvas* canvas) {}
    
};

static void draw_hilite(GCanvas* canvas, Shape* shape) {
    canvas->save();
    shape->concat(canvas);

    GRect r = shape->getRect();
    const float size = CORNER_SIZE;
    GColor c = GColor::MakeARGB(1, 0, 0, 0);
    draw_corner(canvas, c, r.fLeft, r.fTop, size, size);
    draw_corner(canvas, c, r.fLeft, r.fBottom, size, -size);
    draw_corner(canvas, c, r.fRight, r.fTop, -size, size);
    draw_corner(canvas, c, r.fRight, r.fBottom, -size, -size);

    canvas->restore();
}

void Shape::drawHilite(GCanvas* canvas) {
    draw_hilite(canvas, this);
}

static void constrain_color(GColor* c) {
    c->fA = std::max(std::min(c->fA, 1.f), 0.1f);
    c->fR = std::max(std::min(c->fR, 1.f), 0.f);
    c->fG = std::max(std::min(c->fG, 1.f), 0.f);
    c->fB = std::max(std::min(c->fB, 1.f), 0.f);
}

class ShapeClickHandler : public GClick {
public:
    ShapeClickHandler(GPoint loc) : GClick(loc, "ShapeClickHandler") {}
    
    virtual void handleClick(GWindow*) = 0;
};

class RectShape : public Shape {
public:
    RectShape(GColor c) : fColor(c) {
        fRect = GRect::MakeXYWH(0, 0, 0, 0);
    }

    virtual void onDraw(GCanvas* canvas) {
        canvas->fillRect(fRect, fColor);
    }

    virtual GRect getRect() { return fRect; }
    virtual void setRect(const GRect& r) { fRect = r; }
    virtual GColor getColor() { return fColor; }
    virtual void setColor(const GColor& c) { fColor = c; }

private:
    GRect   fRect;
    GColor  fColor;
};

class BitmapShape : public Shape {
public:
    BitmapShape(const GBitmap& bm) : fBM(bm) {
        const int w = std::max(bm.width(), 100);
        const int h = std::max(bm.height(), 100);
        fRect = GRect::MakeXYWH(100, 100, w, h);
    }
    
    virtual void onDraw(GCanvas* canvas) {
        canvas->fillBitmapRect(fBM, fRect);
    }
    
    virtual GRect getRect() { return fRect; }
    virtual void setRect(const GRect& r) { fRect = r; }
    virtual GColor getColor() { return GColor::MakeARGB(1, 0, 0, 0); }
    virtual void setColor(const GColor&) {}
    
private:
    GRect   fRect;
    GBitmap fBM;
};

class GradientShape : public Shape {
public:
    GradientShape(const GColor& c0, const GColor& c1) : fColors{c0, c1} {
        fRect = GRect::MakeXYWH(100, 100, 200, 100);
        fA = 0.5f;
    }
    
    virtual void onDraw(GCanvas* canvas) {
        const GPoint pts[] = {
            GPoint::Make(fRect.left(), fRect.top()),
            GPoint::Make(fRect.right(), fRect.bottom()),
        };
        GColor colors[2]{ fColors[0], fColors[1] };
        colors[0].fA *= fA * 2;
        colors[1].fA *= fA * 2;
        GShader* shader = GShader::FromLinearGradient(pts, colors);
        if (shader) {
            canvas->shadeRect(fRect, shader);
            delete shader;
        }
    }
    
    virtual GRect getRect() { return fRect; }
    virtual void setRect(const GRect& r) { fRect = r; }
    virtual GColor getColor() { return GColor::MakeARGB(fA, 0, 0, 0); }
    virtual void setColor(const GColor& c) {
        fA = c.fA;
    }

private:
    GRect   fRect;
    GColor  fColors[2];
    float   fA;
};

static void make_regular_poly(GPoint pts[], int count, float cx, float cy, float radius) {
    float angle = 0;
    const float deltaAngle = M_PI * 2 / count;
    
    for (int i = 0; i < count; ++i) {
        pts[i].set(cx + cos(angle) * radius, cy + sin(angle) * radius);
        angle += deltaAngle;
    }
}

class RadialShape : public Shape {
public:
    RadialShape(const GColor& c0, const GColor& c1) : fColors{c0, c1} {
        fRect = GRect::MakeXYWH(100, 100, 200, 100);
        fA = 0.5f;

        make_regular_poly(fPoly, GARRAY_COUNT(fPoly), 0, 0, 0.5);
    }
    
    virtual void onDraw(GCanvas* canvas) {
        canvas->save();
        canvas->translate(fRect.x(), fRect.y());
        canvas->scale(fRect.width(), fRect.height());
        canvas->translate(0.5, 0.5);

        GColor colors[2]{ fColors[0], fColors[1] };
        colors[0].fA *= fA * 2;
        colors[1].fA *= fA * 2;
        GShader* shader = GShader::FromRadialGradient(GPoint::Make(0, 0), 1, colors);
        if (shader) {
            canvas->shadeConvexPolygon(fPoly, GARRAY_COUNT(fPoly), shader);
            delete shader;
        }
        
        canvas->restore();
    }
    
    virtual GRect getRect() { return fRect; }
    virtual void setRect(const GRect& r) { fRect = r; }
    virtual GColor getColor() { return GColor::MakeARGB(fA, 0, 0, 0); }
    virtual void setColor(const GColor& c) {
        fA = c.fA;
    }
    
private:
    GRect   fRect;
    GColor  fColors[2];
    float   fA;
    GPoint  fPoly[64];
};

class StrokeShape : public Shape {
public:
    StrokeShape() : fPts{{ 30, 30 }, {200, 60}, {50, 200}} {
        fRect = GRect::MakeXYWH(100, 100, 200, 100);
        fColor = GColor::MakeARGB(1, 1, 0, 0);
        fWidth = 25;
        fMiterLimit = 5;
        fAddCap = false;
        fIsClosed = false;
    }
    
    void onDraw(GCanvas* canvas) override {
        GShader* shader = GShader::FromColor(fColor);
        GCanvas::Stroke stroke = { fWidth, fMiterLimit, fAddCap };
        canvas->strokePolygon(fPts, 3, fIsClosed, stroke, shader);
        delete shader;
    }
    
    void drawHilite(GCanvas* canvas) override {
        const float r = 3;
        GColor c = GColor::MakeARGB(1, 0, 0, 1);
        for (auto p : fPts) {
            canvas->fillRect(GRect::MakeLTRB(p.fX - r, p.fY - r, p.fX + r, p.fY + r), c);
        }
    }
    
    GRect getRect() override { return fRect; }
    void setRect(const GRect& r) override { fRect = r; }
    GColor getColor() override { return fColor; }
    void setColor(const GColor& c) override { fColor = c; }

    class MyHandler : public ShapeClickHandler {
        GPoint* fPtr;
    public:
        MyHandler(GPoint loc, GPoint* ptr) : ShapeClickHandler(loc), fPtr(ptr) {}

        void handleClick(GWindow* wind) {
            *fPtr = this->curr();
            wind->requestDraw();
        }
    };

    GClick* findClickHandler(GPoint loc) override {
        for (int i = 0; i < 3; ++i) {
            if (hit_test(loc.fX, loc.fY, fPts[i].fX, fPts[i].fY)) {
                return new MyHandler(loc, &fPts[i]);
            }
        }
        return nullptr;
    }

    bool handleKey(int c) override {
        switch (c) {
            case 'W': fWidth += 4; return true;
            case 'w': fWidth -= 4; return true;
            case 'c':
            case 'C': fAddCap = !fAddCap; return true;
            case 'm':
            case 'M': fMiterLimit = -fMiterLimit; return true;
            case 'x':
            case 'X': fIsClosed = !fIsClosed; return true;
        }
        return this->INHERITED::handleKey(c);
    }

private:
    GRect   fRect;
    float   fA;
    GPoint  fPts[3];
    GColor  fColor;
    float   fWidth;
    float   fMiterLimit;
    bool    fAddCap;
    bool    fIsClosed;

    typedef Shape INHERITED;
};

static void alloc_bitmap(GBitmap* bm, int w, int h) {
    bm->fWidth = w;
    bm->fHeight = h;
    bm->fRowBytes = bm->fWidth * sizeof(GPixel);
    bm->fPixels = (GPixel*)calloc(bm->height(), bm->rowBytes());
}

static Shape* cons_up_shape(int index) {
    const char* names[] = {
        "apps/spock.png",
        "expected/blend_black.png",
        "expected/circles_blend.png",
        "expected/solid_ramp.png",
        "expected/spocks_zoom.png",
        "expected/blend_white.png",
        "expected/circles_fat.png",
        "expected/spocks_quad.png",
    };
    GBitmap bm;
    if (index < GARRAY_COUNT(names) && bm.readFromFile(names[index])) {
        return new BitmapShape(bm);
    }
    return NULL;
}

class ResizeClick : public GClick {
    GPoint  fAnchor;
public:
    ResizeClick(GPoint loc, GPoint anchor) : GClick(loc, "resize"), fAnchor(anchor) {}

    GRect makeRect() const { return make_from_pts(this->curr(), fAnchor); }
};

class TestWindow : public GWindow {
    std::vector<Shape*> fList;
    Shape* fShape;
    GColor fBGColor;

public:
    TestWindow(int w, int h) : GWindow(w, h) {
        fBGColor = GColor::MakeARGB(1, 1, 1, 1);
        fShape = NULL;
    }

    virtual ~TestWindow() {}
    
protected:
    virtual void onDraw(GCanvas* canvas) {
        canvas->fillRect(GRect::MakeXYWH(0, 0, 10000, 10000), fBGColor);

        for (int i = 0; i < fList.size(); ++i) {
            fList[i]->draw(canvas);
        }
        if (fShape) {
            fShape->drawHilite(canvas);
        }
    }

    virtual bool onKeyPress(const XEvent&, KeySym sym) {
        if (fShape && fShape->handleKey(sym)) {
            this->updateTitle();
            this->requestDraw();
            return true;
        }

        if (sym >= '1' && sym <= '9') {
            fShape = cons_up_shape(sym - '1');
            if (fShape) {
                fList.push_back(fShape);
                this->updateTitle();
                this->requestDraw();
            }
        }
        if ('d' == sym) {
            fShape = new GradientShape(rand_color(), rand_color());
            fList.push_back(fShape);
            this->updateTitle();
            this->requestDraw();
            return true;
        }
        
        if ('c' == sym) {
            fShape = new RadialShape(rand_color(), rand_color());
            fList.push_back(fShape);
            this->updateTitle();
            this->requestDraw();
            return true;
        }
        if ('s' == sym) {
            fShape = new StrokeShape;
            fList.push_back(fShape);
            this->updateTitle();
            this->requestDraw();
            return true;
        }
        
        if (fShape) {
            switch (sym) {
                case XK_Left:
                    fShape->fAngle -= 5;
                    this->requestDraw();
                    break;

                case XK_Right:
                    fShape->fAngle += 5;
                    this->requestDraw();
                    break;

                case XK_Up: {
                    int index = find_index(fList, fShape);
                    if (index < fList.size() - 1) {
                        std::swap(fList[index], fList[index + 1]);
                        this->requestDraw();
                        return true;
                    }
                    return false;
                }
                case XK_Down: {
                    int index = find_index(fList, fShape);
                    if (index > 0) {
                        std::swap(fList[index], fList[index - 1]);
                        this->requestDraw();
                        return true;
                    }
                    return false;
                }
                case XK_BackSpace:
                    this->removeShape(fShape);
                    fShape = NULL;
                    this->updateTitle();
                    this->requestDraw();
                    return true;
                default:
                    break;
            }
        }

        GColor c = fShape ? fShape->getColor() : fBGColor;
        const float delta = 0.1f;
        switch (sym) {
            case 'a': c.fA -= delta; break;
            case 'A': c.fA += delta; break;
            case 'r': c.fR -= delta; break;
            case 'R': c.fR += delta; break;
            case 'g': c.fG -= delta; break;
            case 'G': c.fG += delta; break;
            case 'b': c.fB -= delta; break;
            case 'B': c.fB += delta; break;
            default:
                return false;
        }
        constrain_color(&c);
        if (fShape) {
            fShape->setColor(c);
        } else {
            c.fA = 1;   // need the bg to stay opaque
            fBGColor = c;
        }
        this->updateTitle();
        this->requestDraw();
        return true;
    }

    virtual GClick* onFindClickHandler(GPoint loc) {
        if (fShape) {
            GClick* handler = fShape->findClickHandler(loc);
            if (handler) {
                return handler;
            }

            GPoint anchor;
            if (in_resize_corner(fShape->getRect(), loc.x(), loc.y(), &anchor)) {
                return new ResizeClick(loc, anchor);
            }
        }

        for (int i = fList.size() - 1; i >= 0; --i) {
            if (contains(fList[i]->getRect(), loc.x(), loc.y())) {
                fShape = fList[i];
                this->updateTitle();
                return new GClick(loc, "move");
            }
        }
        
        // else create a new shape
        fShape = new RectShape(rand_color());
        fList.push_back(fShape);
        this->updateTitle();
        return new GClick(loc, "create");
    }

    virtual void onHandleClick(GClick* click) {
        if (click->isName("ShapeClickHandler")) {
            ((ShapeClickHandler*)click)->handleClick(this);
            return;
        }
        
        if (click->isName("move")) {
            const GPoint curr = click->curr();
            const GPoint prev = click->prev();
            fShape->setRect(offset(fShape->getRect(), curr.x() - prev.x(), curr.y() - prev.y()));
        } else if (click->isName("resize")) {
            fShape->setRect(((ResizeClick*)click)->makeRect());
        } else {
            fShape->setRect(make_from_pts(click->orig(), click->curr()));
        }
        if (NULL != fShape && GClick::kUp_State == click->state()) {
            if (fShape->getRect().isEmpty()) {
                this->removeShape(fShape);
                fShape = NULL;
            }
        }
        this->updateTitle();
        this->requestDraw();
    }

private:
    void removeShape(Shape* target) {
        GASSERT(target);

        std::vector<Shape*>::iterator it = std::find(fList.begin(), fList.end(), target);
        if (it != fList.end()) {
            fList.erase(it);
        } else {
            GASSERT(!"shape not found?");
        }
    }

    void updateTitle() {
        char buffer[100];
        buffer[0] = ' ';
        buffer[1] = 0;

        GColor c = fBGColor;
        if (fShape) {
            c = fShape->getColor();
        }

        sprintf(buffer, "%02X %02X %02X %02X",
                int(c.fA * 255), int(c.fR * 255), int(c.fG * 255), int(c.fB * 255));
        this->setTitle(buffer);
    }

    typedef GWindow INHERITED;
};

int main(int argc, char const* const* argv) {
    GWindow* wind = new TestWindow(640, 480);

    return wind->run();
}

