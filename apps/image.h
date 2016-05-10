/**
 *  Copyright 2015 Mike Reed
 *
 *  575 Introduction to Graphics
 */

#ifndef CS575_image_DEFINED
#define CS575_image_DEFINED

class GCanvas;

struct CS575DrawRec {
    void        (*fDraw)(GCanvas*);
    int         fWidth;
    int         fHeight;
    int         fPANumber;
    const char* fName;
};

/*
 *  Array is terminated when fDraw is NULL
 */
extern const CS575DrawRec gDrawRecs[];

#endif
