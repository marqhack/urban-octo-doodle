/**
 *  Copyright 2015 Mike Reed
 */

#ifndef GBitmap_DEFINED
#define GBitmap_DEFINED

#include "GPixel.h"

class GBitmap {
public:
    int width() const { return fWidth; }
    int height() const { return fHeight; }
    size_t rowBytes() const { return fRowBytes; }
    GPixel* pixels() const { return fPixels; }

    void reset() {
        fWidth = 0;
        fHeight = 0;
        fPixels = NULL;
        fRowBytes = 0;
    }

    int     fWidth;
    int     fHeight;
    GPixel* fPixels;
    size_t  fRowBytes;

    GPixel* getAddr(int x, int y) const {
        GASSERT(x >= 0 && x < this->width());
        GASSERT(y >= 0 && y < this->height());
        return this->pixels() + x + (y * this->rowBytes() >> 2);
    }

    /**
     *  Attempt to read the png image stored in the named file.
     *
     *  On success, allocate the memory for the pixels using malloc() and set bitmap to the result,
     *  returning true. The caller must call free(bitmap->fPixels) when they are finished.
     *
     *  On failure, return false and bitmap is reset to empty.
     */
    bool readFromFile(const char path[]);

    /*
     *  Attempt to write the bitmap as a PNG into a new file (the file will be created/overwritten).
     *  Return true on success.
     */
    bool writeToFile(const char path[]) const;
};

#endif
