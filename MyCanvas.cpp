#include "GCanvas.h"
#include "GShader.h"
#include "GBitmap.h"
#include "GColor.h"
#include "GRect.h"
#include <iostream>
#include "GMath.h"
#include "GPixel.h"
#include "GPoint.h"
#include <cmath> 
#include <stack>

using namespace std;

struct Edge {
  GPoint A; 
  GPoint B;
};
float** createIdentityMatrix(int rows, int columns);     
float** multiplyMatricies(float** A, float** B);
Edge* getBoundaryEdge(Edge* edges [], int num_edges, Edge* current_edge);
int computeEdgeBoundary(Edge* edge, int y);
void clipBottom(Edge* edges [], int i, GRect rect);
void clipTop(Edge* edges [], int i, GRect rect);
void clipLeft(Edge* edges [], int i, GRect rect, int count);
void clipRight(Edge* edges [], int i, GRect rect, int count);
bool isSamePoint(GPoint A, GPoint B);
GPoint getLowerPoint(Edge* edge);
GPoint getUpperPoint(Edge* edge);
GPoint getRightPoint(Edge* edge);
GPoint getLeftPoint(Edge* edge);
void sortEdges(Edge* edges [], int count);
void printEdges(Edge* edges [], int count);
void printPts(const GPoint pts [], int count);
void printMatrix(float** matrix);
GPixel lookupPixel(float x, float y, float** transform, const GBitmap& bitmap);
float** createTransformMatrix(const GBitmap& bitmap, const GRect& rect);
int clamp(int min, int value, int max);
float** createSrcDst(const GBitmap& bitmap, const GRect& rect);
float** createScaleMatrix(const GBitmap& bitmap, const GRect rect);
float** createTranslateMatrix(const GBitmap& bitmap, const GRect rect);
float getDstY(float x, float y, float** transform);
float getDstX(float x, float y, float** transform);
GPixel srcOverBlend(GPixel src, GPixel dst);
GPixel makePixelFromColor(GColor color);
float** invertMatrix(float** matrix);
float getUnitValue(int x, int y, float** M);
float** mapToLine(const GPoint* pts);
GPoint projectPointontoLine(const GPoint* pts, GPoint kpt);
GColor linearCombineColors(float unitValue, const GColor* colors);
float distanceBetween(GPoint A, GPoint B);
float getVectorLength(GPoint A);
GPoint translatePoint(GPoint pt, float** matrix);
GPoint makeUnitVector(GPoint v);
float getMiterLength(GPoint u, GPoint v, float width, float limit);
GPoint vectorFromEdge(Edge* e);

/* The constructor for creating a my canvas object, subclassed from GCanvas 
 * Implements the virtual methods in GCanvas 
 * associates a GRect with the bitmap;
*/
class MyCanvas : public GCanvas{
public:

  MyCanvas(const GBitmap& bitmap){
    the_bitmap = bitmap;
    the_rect = GRect::MakeLTRB((float)0, (float)0, (float)the_bitmap.fWidth, (float)the_bitmap.fHeight);
    ctm = createIdentityMatrix(3,3);
    save();
  }

  void clear(const GColor& color);
  void fillRect(const GRect& rect, const GColor& color);
  void fillBitmapRect(const GBitmap& src , const GRect& dst);
  void fillConvexPolygon(const GPoint pts [], int count, const GColor&);
  void fillConvexPolygonBitmap(Edge* edges [], int count, const GBitmap& src, float** untransform);
  void save();
  void restore();
  void concat(const float matrix[6]);
  void shadeRect(const GRect& rect, GShader* shader);
  void shadeConvexPolygon(const GPoint[], int count, GShader* shader);
  void strokePolygon(const GPoint[], int n, bool isClosed, const Stroke&, GShader*);

protected:
  stack <float**> ctmStack;
  GBitmap the_bitmap;
  GRect the_rect;
  float** ctmSave = createIdentityMatrix(3, 3);
  float** ctm = createIdentityMatrix(3, 3);

};


/* This is the public facing function to create a new GCanvas object
 *   It checks to make sure that the bitmap passed in is all kosher.
 *   Calls MyCanvas constructor, which adds a GRect to simply things like
       intersection, etc.
*/
GCanvas* GCanvas::Create(const GBitmap& bitmap){
  if(bitmap.fHeight <= 0){return NULL;} 
  if(bitmap.fWidth <= 0){return NULL;}
  if(bitmap.fRowBytes < bitmap.fWidth * sizeof(GPixel)){return NULL;}
  if(bitmap.fPixels == NULL){return NULL;}
  return new MyCanvas(bitmap); 
}


/******************************************************************************
                  Helper Methods
******************************************************************************/
/* Creates identity matrix of size rows x columns */
float** createIdentityMatrix(int rows, int columns){     
    float** matrix = 0;    
    // Make identity matrix 
    matrix = new float*[rows];
    for(int r = 0; r < rows; r++){
        matrix[r] = new float[columns];
        for(int c = 0; c < columns; c++){
            if(c == r){ matrix[r][c] = 1; }
            else{ matrix[r][c] = 0; }
        } 
    }
    return matrix;
}
/* Takes in two pixels, does Porter-Duffer src-over mode blending, and returns a blended pixel*/
GPixel srcOverBlend(GPixel src, GPixel dst){

  // Get components of src color 
  unsigned int srcA = GPixel_GetA(src);       
  unsigned int srcR = GPixel_GetR(src);       
  unsigned int srcG = GPixel_GetG(src);       
  unsigned int srcB = GPixel_GetB(src);       
 
  // Get components of dst color 
  unsigned int dstA =  GPixel_GetA(dst);       
  unsigned int dstR =  GPixel_GetR(dst);       
  unsigned int dstG =  GPixel_GetG(dst);       
  unsigned int dstB =  GPixel_GetB(dst);       
  
  /* Blend using src over, Do the multiply thing */
  unsigned int tA = srcA + (dstA - ( ( ((srcA*dstA)*65793+(1<<23))>>24 )   ));
  unsigned int tR = srcR + ( ( dstR*(255-srcA) *65793+(1<<23) )>>24 );
  unsigned int tG = srcG + ( ( dstG*(255-srcA) *65793+(1<<23) )>>24 );
  unsigned int tB = srcB + ( ( dstB*(255-srcA) *65793+(1<<23) )>>24 );
        
  GPixel pixel = GPixel_PackARGB(tA, tR, tG, tB);
  return pixel;
}
/* Map rect coordinates onto bitmap coordinates using matrix
   Analytical formula instead of matrix multiplications
*/
float** createTransformMatrix(const GBitmap& bitmap, const GRect& rect){
 
  float** transform = createIdentityMatrix(3, 3);
  
  float Sx = bitmap.width()/rect.width();
  float Sy = bitmap.height()/rect.height();

  float Tx = 0-rect.left();
  float Ty = 0-rect.top();

  transform[0][0] = Sx;
  transform[1][1] = Sy;

  transform[0][2] = Sx*Tx; 
  transform[1][2] = Sy*Ty; 


 
  return transform;
}

/* takes in a coordinate (x, y), and a transform matrix, and it finds the corresponding transformed pixel from bitmap, and returns it */ 
GPixel lookupPixel(float x, float y, float** transform, const GBitmap& bitmap){

  int x_coord = x*transform[0][0] + y*transform[0][1] + transform[0][2];
  int y_coord = x*transform[1][0] + y*transform[1][1] + transform[1][2];
  int z_coord = x*transform[2][0] + y*transform[2][1] + transform[2][2];

  int lookup_y = clamp(0, (int)y_coord, bitmap.height() - 1);
  int lookup_x = clamp(0, (int)x_coord, bitmap.width() -1);


  return bitmap.getAddr(lookup_x, lookup_y)[0]; 
 

}


/******************************************************************************
                            PA1
*******************************************************************************/
void MyCanvas::clear(const GColor& color){
  /* Get values from color */
  float fA = color.fA;
  float fR = color.fR;
  float fG = color.fG;
  float fB = color.fB;     
  
  /* Pre-multiply values by alpha */
  fR = fA*fR;
  fG = fA*fG;
  fB = fA*fB; 
    
  /* scale to 255 */
  unsigned int rA = fA *255; 
  unsigned int rR = fR *255;
  unsigned int rG = fG *255; 
  unsigned int rB = fB *255;
  
  
  GPixel* dst = the_bitmap.fPixels;
  for (int y = 0; y < the_bitmap.height(); ++y) {
    for (int x = 0; x < the_bitmap.width(); ++x) {

      dst[x] = GPixel_PackARGB(rA, rR, rG, rB); 

    }
    dst = (GPixel*)((char*)dst + the_bitmap.rowBytes());
  }

}

/* fills a rectangle using porter-duff src-over mode */
void MyCanvas::fillRect(const GRect& rect, const GColor& color)
{
   
  /* Find intersection of rectangle */
  GRect o_rect = GRect::MakeLTRB(rect.fLeft, rect.fTop, rect.fRight, rect.fBottom);
  o_rect.intersect(the_rect);
  GIRect other_rect = o_rect.round(); // this is the rectange we draw into.
  
  /* Get components of src color */
  float fA = color.fA;
  float fR = color.fR;
  float fG = color.fG;
  float fB = color.fB;    
  if(fA > 1 || fA < 0){ return;} 
  if(fR > 1 || fR < 0){ return;} 
  if(fG > 1 || fG < 0){ return;} 
  if(fB > 1 || fB < 0){ return;} 


  /* Pre-multiply values by alpha */
  /* scale to 255 (These are the premulled components of src pixel) */
  unsigned int srcA = fA *255; 
  unsigned int srcR = fA*fR *255;
  unsigned int srcG = fA*fG *255; 
  unsigned int srcB = fA*fB *255;
  
  /* pack a pixel with this color, because we have a method that blends two pixel */
  GPixel src_pixel =  GPixel_PackARGB(srcA, srcR, srcG, srcB);

  
  GPixel* dst = the_bitmap.fPixels+other_rect.fLeft;
  dst = (GPixel*)((char*)dst + the_bitmap.rowBytes()*other_rect.fTop); 
  for (int y = 0; y < other_rect.height(); ++y) {
    for (int x = 0; x < other_rect.width(); ++x) {
      GPixel dst_pixel = dst[x];
      dst[x] = srcOverBlend(src_pixel, dst_pixel);
    }
    dst = (GPixel*)((char*)dst + the_bitmap.rowBytes());
  }
}
 

/******************************************************************************
                            PA2
*******************************************************************************/
/* scale a bitmap (src) to fill the rectangle described by dst
*/
void MyCanvas::fillBitmapRect(const GBitmap& src , const GRect& dst){
     
  float** src_to_dst = createSrcDst(src,dst); 


  float** untransform = invertMatrix(multiplyMatricies(ctm,src_to_dst));
  
  float** id = createIdentityMatrix(3,3); 
  GRect other_rect = dst; 

  int x0 =  getDstX(other_rect.left(), other_rect.top(),this->ctm); 
  int y0 =  getDstY(other_rect.left(), other_rect.top(),this->ctm); 

  int x1 =  getDstX(other_rect.right(), other_rect.top(),this->ctm); 
  int y1 =  getDstY(other_rect.right(), other_rect.top(),this->ctm); 

  int x2 =  getDstX(other_rect.right(), other_rect.bottom(),this->ctm); 
  int y2 =  getDstY(other_rect.right(), other_rect.bottom(),this->ctm); 

  int x3 =  getDstX(other_rect.left(), other_rect.bottom(),this->ctm); 
  int y3 =  getDstY(other_rect.left(), other_rect.bottom(),this->ctm); 

  GPoint pts [4];
  pts[3] = GPoint::Make(x0,y0);
  pts[2] = GPoint::Make(x1,y1);
  pts[1] = GPoint::Make(x2,y2);
  pts[0] = GPoint::Make(x3,y3);

  const int count = 4;
  Edge* edges [count*3];

  for(int i = 0; i < count*3; i++){
    edges[i] = NULL;
  } 

  for(int i = 0; i < count; i++){
    Edge* edge = new Edge;
    edge->A = pts[i];
    edge->B = pts[(i+1) % count]; 
    edges[i] = edge;
  }

  for(int i = 0; i < count; i++){
    clipTop(edges, i, the_rect); 
    clipBottom(edges, i, the_rect);
    clipLeft(edges, i, the_rect, count);
    clipRight(edges, i, the_rect, count); 
  }

  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(edges[i]->A.y() < the_rect.top() && edges[i]->B.y() < the_rect.top()){
        edges[i] = NULL;
      }else if(edges[i]->A.y() > the_rect.bottom() && edges[i]->B.y() > the_rect.bottom()){
        edges[i] = NULL;
      }
    }
  }

  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(edges[i]->A.y() == edges[i]->B.y()){
        edges[i] = NULL;

      }
    }
  }

  fillConvexPolygonBitmap(edges, 4, src, untransform);
}


/******************************************************************************
                 PA3 
******************************************************************************/

/**
 *  Fill the convex polygon with the color, following the same "containment" rule as
 *  rectangles.
 *
 *  Any area in the polygon that is outside of the bounds of the canvas is ignored.
 *
 *  If the color's alpha is < 1, blend it using SRCOVER blend mode.
 */
void MyCanvas::fillConvexPolygon(const GPoint pts [], int count, const GColor& color){

  if(count < 1){return;}

  GPoint* ctmPts = new GPoint[count];
 
  
  for(int i = 0; i < count; i++){     
    int x_trans;
    int y_trans;
    x_trans = getDstX(pts[i].x(), pts[i].y(),this->ctm); 
    y_trans = getDstY(pts[i].x(), pts[i].y(),this->ctm);
    ctmPts[i] = GPoint::Make(x_trans, y_trans);
  }

  Edge** edges = new Edge*[count*3];

  //fill array with edges/NULL spaces
  for(int i = 0; i < count*3; i++){
    edges[i] = NULL;
  } 
  // create initial, unclipped polygon edges
  for(int i = 0; i < count; i++){
    Edge* edge = new Edge;
    edge->A = ctmPts[i];
    edge->B = ctmPts[(i+1) % count]; 
    edges[i] = edge;
  }


  //Clip edges 
  for(int i = 0; i < count; i++){
    clipTop(edges, i, the_rect); 
    clipBottom(edges, i, the_rect);
    clipLeft(edges, i, the_rect, count);
    clipRight(edges, i, the_rect, count);
  }

  // discard edges completely above top, or below bottom 
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(edges[i]->A.y() <= the_rect.top() && edges[i]->B.y() <= the_rect.top()){
        edges[i] = NULL;
      }else if(edges[i]->A.y() >= the_rect.bottom() && edges[i]->B.y() >= the_rect.bottom()){
        edges[i] = NULL;
      }
    }
  }


  // get rid of horizontal edges
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(edges[i]->A.y() == edges[i]->B.y()){
        edges[i] = NULL;

      }
    }
  }

  //Find highest and lowest points
  GPoint maxYPt = GPoint::Make(the_rect.left(), the_rect.bottom()); // closest to top
  GPoint minYPt = GPoint::Make(the_rect.left(), the_rect.top()); // closest to bottom
 
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(getUpperPoint(edges[i]).y() < maxYPt.y()){ 
        maxYPt = getUpperPoint(edges[i]);
      }
      if(getLowerPoint(edges[i]).y() > minYPt.y()){
        minYPt = getLowerPoint(edges[i]);
      }
    }
  }
 
  // these will be the two highest hedges
  Edge* first_edge = new Edge;
  Edge* second_edge = new Edge;

  int highestEdgeIndex = -1;
  // get highest edge
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(getUpperPoint(edges[i]).y() <= maxYPt.y()){  
        first_edge = edges[i];
        highestEdgeIndex = i;
      }
    }
  }

  //second edge = minedge
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(getLowerPoint(edges[i]).y() >= minYPt.y()){  
        second_edge = edges[i];
      }
    }
  }
  
  // get second highest edge
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(getUpperPoint(edges[i]).y() <= getUpperPoint(second_edge).y()){        
        if(i != highestEdgeIndex){  
          second_edge = edges[i];
        }
      }
    }
  }

  //see which is left edge, which is right edge
  Edge* left_edge = new Edge;
  Edge* right_edge = new Edge;

  if(getLowerPoint(first_edge).x() < getLowerPoint(second_edge).x()){
    left_edge = first_edge;
    right_edge = second_edge;
  }else{
    left_edge = second_edge;
    right_edge = first_edge;    
  }

  for(int y = maxYPt.y(); y < minYPt.y(); y++){
  
    if( y > getLowerPoint(left_edge).y()){ 
      left_edge = getBoundaryEdge(edges, count, left_edge);
    }
    if( y > getLowerPoint(right_edge).y()){
      right_edge = getBoundaryEdge(edges, count, right_edge);
    } 
    int leftBoundary = computeEdgeBoundary(left_edge, y);
    int rightBoundary = computeEdgeBoundary(right_edge, y);
  
    GPixel* dst = the_bitmap.fPixels;
    dst = (GPixel*)((char*)dst + the_bitmap.rowBytes()*y); 

    for(int x = leftBoundary; x < rightBoundary; x++){
      GPixel src_pixel = makePixelFromColor(color); 
      GPixel dst_pixel = dst[x];
      dst[x] = srcOverBlend(src_pixel, dst_pixel);
    }  
    
  }

}

/******************************************************************************
                 PA4 
******************************************************************************/
/**
 *  Saves a copy of the CTM, allowing subsequent modifications (by calling concat()) to be
 *  undone when restore() is called.
 *
*/
void MyCanvas::save(){

  this->ctmStack.push(this->ctm);

/*
  float** m = createIdentityMatrix(3,3);

  m[0][0] = this->ctm[0][0];
  m[0][1] = this->ctm[0][1];
  m[0][2] = this->ctm[0][2];

  m[1][0] = this->ctm[1][0];
  m[1][1] = this->ctm[1][1];
  m[1][2] = this->ctm[1][2];

  this->ctmSave = m;
*/

}

/**
 *  Balances calls to save(), returning the CTM to the state it was in when the corresponding
 *  call to save() was made. These calls can be nested.
 *
*/
void MyCanvas::restore(){

  this->ctm = this->ctmStack.top(); 
  this->ctmStack.pop(); 
/*
  float** m = createIdentityMatrix(3,3);

  m[0][0] = this->ctmSave[0][0];
  m[0][1] = this->ctmSave[0][1];
  m[0][2] = this->ctmSave[0][2];

  m[1][0] = this->ctmSave[1][0];
  m[1][1] = this->ctmSave[1][1];
  m[1][2] = this->ctmSave[1][2];
  
  this->ctm = m; 
*/

}

/**
 *  Modifies the CTM (current transformation matrix) by pre-concatenating it with the specfied
 *  matrix.
 *
 *  CTM' = CTM * matrix
*/
void MyCanvas::concat(const float matrix[6]){
 
  float** m = createIdentityMatrix(3,3);
  // because its only 6 elements, we use the analytical formula 
  m[0][0] = matrix[0];
  m[0][1] = matrix[1];
  m[0][2] = matrix[2];

  m[1][0] = matrix[3];
  m[1][1] = matrix[4];
  m[1][2] = matrix[5];

  this->ctm = multiplyMatricies(ctm, m);

}

/******************************************************************************
                 PA5 
******************************************************************************/
void MyCanvas::shadeRect(const GRect& rect, GShader* shader){

  GPoint pts [4];
  pts[3] = GPoint::Make(rect.left(), rect.top());
  pts[2] = GPoint::Make(rect.right(), rect.top());
  pts[1] = GPoint::Make(rect.right(), rect.bottom());
  pts[0] = GPoint::Make(rect.left(), rect.bottom());

  shadeConvexPolygon(pts, 4, shader); 
  
}

void MyCanvas::shadeConvexPolygon(const GPoint pts [], int count, GShader* shader){
  GPoint* ctmPts = new GPoint[count];

  float x_trans;
  float y_trans;
  
  for(int i = 0; i < count; i++){ 
    x_trans = getDstX(pts[i].x(), pts[i].y(),this->ctm); 
    y_trans = getDstY(pts[i].x(), pts[i].y(),this->ctm);
    ctmPts[i] = GPoint::Make(x_trans, y_trans);
  }
  
  Edge** edges = new Edge*[count*3];
  for(int i = 0; i < count*3; i++){
    edges[i] = NULL;
  } 
  for(int i = 0; i < count; i++){
    Edge* edge = new Edge;
    edge->A = ctmPts[i];
    edge->B = ctmPts[(i+1) % count]; 
    edges[i] = edge;
  } 
   
  for(int i = 0; i < count; i++){
    clipTop(edges, i, the_rect); 
    clipBottom(edges, i, the_rect);
    clipLeft(edges, i, the_rect, count);
    clipRight(edges, i, the_rect, count);
  }
  // discard edges completely above top, or below bottom 
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(edges[i]->A.y() <= the_rect.top() && edges[i]->B.y() <= the_rect.top()){
        edges[i] = NULL;
      }else if(edges[i]->A.y() >= the_rect.bottom() && edges[i]->B.y() >= the_rect.bottom()){
        edges[i] = NULL;
      }
    }
  }
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(edges[i]->A.y() == edges[i]->B.y()){
        edges[i] = NULL;
      }
    }
  }

  
 
  


  GPoint maxYPt = GPoint::Make(the_rect.left(), the_rect.bottom()); // closest to top
  GPoint minYPt = GPoint::Make(the_rect.left(), the_rect.top()); // closest to bottom 
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(getUpperPoint(edges[i]).y() < maxYPt.y()){ 
        maxYPt = getUpperPoint(edges[i]);
      }
      if(getLowerPoint(edges[i]).y() > minYPt.y()){
        minYPt = getLowerPoint(edges[i]);
      }
    }
  }
 
  // these will be the two highest hedges
  Edge* first_edge = new Edge;
  Edge* second_edge = new Edge;

  int highestEdgeIndex = -1;
  // get highest edge
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(getUpperPoint(edges[i]).y() <= maxYPt.y()){  
        first_edge = edges[i];
        highestEdgeIndex = i;
      }
    }
  }

  //second edge = minedge
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(getLowerPoint(edges[i]).y() >= minYPt.y()){  
        second_edge = edges[i];
      }
    }
  }
  
  // get second highest edge
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(getUpperPoint(edges[i]).y() <= getUpperPoint(second_edge).y()){        
        if(i != highestEdgeIndex){  
          second_edge = edges[i];
        }
      }
    }
  }

  //see which is left edge, which is right edge
  Edge* left_edge = new Edge;
  Edge* right_edge = new Edge;

  if(getLowerPoint(first_edge).x() < getLowerPoint(second_edge).x()){
    left_edge = first_edge;
    right_edge = second_edge;
  }else{
    left_edge = second_edge;
    right_edge = first_edge;    
  }

  const float reed[6] = {this->ctm[0][0], this->ctm[0][1], this->ctm[0][2],
  this->ctm[1][0], this->ctm[1][1], this->ctm[1][2]};
  shader->setContext(reed);

  for(int y = maxYPt.y(); y < minYPt.y(); y++){
  
    if( y > getLowerPoint(left_edge).y()){ 
      left_edge = getBoundaryEdge(edges, count, left_edge);
    }
    if( y > getLowerPoint(right_edge).y()){ 
      right_edge = getBoundaryEdge(edges, count, right_edge);
    } 
    int leftBoundary = computeEdgeBoundary(left_edge, y+.5);
    int rightBoundary = computeEdgeBoundary(right_edge, y+.5);
    
    
    GPixel* dst = the_bitmap.fPixels + leftBoundary;
    dst = (GPixel*)((char*)dst + the_bitmap.rowBytes()*y); 

    
    shader->shadeRow(leftBoundary, y, rightBoundary-leftBoundary+1,dst);
   
    



  }

  const float clr[6] = {0,0,1,0,1,0};

  shader->setContext(clr);
}
                 
void clipTop(Edge* edges [], int i, GRect rect){

  float top = rect.top();
  Edge* edge = edges[i];
  

  // nullify after clipping
  if(edges[i]->A.y() <= rect.top() && edges[i]->B.y() <= rect.top()){
    return;
  }
  // don't clip edges already within rect
  if(edge->A.y() >= top && edge->B.y() >= top){
    return;
  }
   
  GPoint p0 = getUpperPoint(edge);
  GPoint p1 = getLowerPoint(edge);

  float t1_height = p1.y() - p0.y();
  float t2_height = top - p0.y();
  float t1_base = p1.x() - p0.x();
  float t2_base = t1_base * (t2_height/t1_height);

  GPoint newPoint = GPoint::Make(p0.x() + t2_base, top);
  
  p0 = newPoint;

  edges[i] = new Edge;
  edges[i]->A = p0;
  edges[i]->B = p1;
}


/* Returns an edge that is within the bottom of the rect 
*/
void clipBottom(Edge* edges [], int i, GRect rect){
  float bottom = rect.bottom();
  Edge* edge = edges[i];

  // edge completely outside: return (discard later)
  if(edges[i]->A.y() >= rect.bottom() && edges[i]->B.y() >= rect.bottom()){
      return;
  }

  // edge completely inside: don't clip
  if(edge->A.y() <= bottom && edge->B.y() <= bottom){ 
    return;
  }

  GPoint p0 = getUpperPoint(edge); 
  GPoint p1 = getLowerPoint(edge);

  float t1_height = p1.y() - p0.y();
  float t2_height = bottom - p0.y();
  float t1_base = p1.x() - p0.x();
  float t2_base = t1_base * (t2_height/t1_height);

  GPoint newPoint = GPoint::Make(p0.x() + t2_base, bottom);
  p1 = newPoint;

  edges[i] = new Edge;
  edges[i]->A = p0;
  edges[i]->B = p1;
  return; 
}



//clips edge (changes contents of edges) and returns any possible new edge
void clipLeft(Edge* edges[], int i, GRect rect, int count){

  float left = rect.left();
  Edge* edge = edges[i];
  

  //edges completely left of left: project onto side
  if(edge->A.x() <= left && edge->B.x() <= left){
    Edge* newEdge = new Edge;
    newEdge->A = GPoint::Make(left, edge->A.y());
    newEdge->B = GPoint::Make(left, edge->B.y());
    delete edges[i];
    edges[i] = newEdge;
    return; 
  }

  // edge completely inside: don't clip
  if(edge->A.x() > left && edge->B.x() > left){
    return;
  }

  GPoint p0 = getLeftPoint(edge);
  GPoint p1 = getRightPoint(edge);

  float t1_height = p1.x() - p0.x();
  float t2_height = left - p0.x();
  float t1_base = p1.y() - p0.y();
  float t2_base = t1_base * (t2_height/t1_height);

  GPoint newPoint = GPoint::Make(left, p0.y()+t2_base);
  Edge* newEdge = new Edge;
  newEdge->A = GPoint::Make(newPoint.x(), p0.y());
  newEdge->B = newPoint;

  p0 = newPoint;
  
  edges[i]->A = p0;
  edges[i]->B = p1;
  edges[i+count] = newEdge;
}

//clips edge (changes contents of edges) and returns any possible new edge
void clipRight(Edge* edges[], int i, GRect rect, int count){

  float right = rect.right();
  Edge* edge = edges[i];

  // edge completely right of right: project onto right
  if(edge->A.x() >= right && edge->B.x() >= right){ 
    Edge* newEdge = new Edge;
    newEdge->A = GPoint::Make(right, edge->A.y());
    newEdge->B = GPoint::Make(right, edge->B.y());
    delete edges[i];
    edges[i] = newEdge;
    return;
  }

  // edge completely inside: don't clip
  if(edge->A.x() < right && edge->B.x() < right){
    return;
  }


  GPoint p0 = getLeftPoint(edges[i]); 
  GPoint p1 = getRightPoint(edges[i]);



  float t1_height = p1.x() - p0.x();
  float t2_height = right - p0.x();
  float t1_base = p1.y() - p0.y();
  float t2_base = t1_base * (t2_height/t1_height);

  GPoint newPoint = GPoint::Make(right, p0.y()+t2_base);
  Edge* newEdge = new Edge;
  newEdge->A = GPoint::Make(newPoint.x(), p1.y());
  newEdge->B = newPoint;

  p1 = newPoint;
 
  edges[i]->A = p0;
  edges[i]->B = p1;
  edges[i+(2*count)] = newEdge;
}

/* Given a current boundary edge, returns a new boundary edge
*/
Edge* getBoundaryEdge(Edge* edges [], int num_edges, Edge* current_edge){ 
  //We always draw from higher y point to lower y point, so our "endpoint" is lower y point
  GPoint startpoint = getUpperPoint(current_edge);
  GPoint endpoint = getLowerPoint(current_edge);
  
  for(int i = 0; i < num_edges*3; i ++){
    if(edges[i] != NULL){
      if( isSamePoint(getUpperPoint(edges[i]), endpoint)){
        
        
        return edges[i];
      }
    } 
  }      
}

/* Returns true if A and B describe the same point
*/
bool isSamePoint(GPoint A, GPoint B){
  return ( ( A.x() == B.x() ) && ( A.y() == B.y() ) );
}

/* Given an edge, and a y coordinate it computes a line equation for the edge to find the corresponding x coord
*/
int computeEdgeBoundary(Edge* edge, int y){
   
  // m = (y1-y0)/(x1-x0)  
  // y = mx + b
  // b = y - mx
  // x = (y-b)/m 

  // if its a vertical line, return x;
  if(edge->A.x() == edge->B.x()) {return edge->A.x();}

  float m = ( edge->A.y() - edge->B.y() ) / (edge->A.x() - edge->B.x());
  float b = edge->A.y() - m*edge->A.x(); 
  float x= ( y-b)/m;

  return floor(x+.5);
}


/* takes a GColor, makes a pixel */
GPixel makePixelFromColor(GColor color){
  // Get values from color 
  float fA = color.fA;
  float fR = color.fR;
  float fG = color.fG;
  float fB = color.fB;     
  
  // Pre-multiply values by alpha 
  fR = fA*fR;
  fG = fA*fG;
  fB = fA*fB; 
   
  // scale to 255 
  unsigned int rA = fA *255; 
  unsigned int rR = fR *255;
  unsigned int rG = fG *255; 
  unsigned int rB = fB *255;
 
  return GPixel_PackARGB(rA, rR, rG, rB); 
}

/* prints a list of edges, for debugging */
void printEdges(Edge* edges [], int count){

  cout << "*********************\n";
  for(int i = 0; i< count*3; i ++){   
    if(edges[i] != NULL){    
    cout << edges[i]->A.x();
    cout << " ";
    cout << edges[i]->A.y();
    cout << "\n";
    cout << edges[i]->B.x();
    cout << " ";
    cout << edges[i]->B.y();
    cout << "\n";
    cout << "\n";
    }
  }
  cout << "\n";
  cout << "\n";
  cout << "\n";
}

void printPts(const GPoint pts [], int count){
  cout << "#############################\n";
  for(int i = 0; i< count; i ++){
   
       
    cout << pts[i].x();
    cout << " ";
    cout << pts[i].y();
    cout << "\n";
  }
  cout << "\n";
}

/* returns point closest to bottom */
GPoint getLowerPoint(Edge* edge){

  if(edge->A.y() > edge->B.y()){
    return edge->A;
  }
  return edge->B;

}

/* returns point closest to top */
GPoint getUpperPoint(Edge* edge){

  if(edge->A.y() < edge->B.y()){
    return edge->A;
  }
  return edge->B;
}

GPoint getLeftPoint(Edge* edge){
  if(edge->A.x() < edge->B.x()){
    return edge->A;
  }
  return edge->B;
}

GPoint getRightPoint(Edge* edge){
  if(edge->A.x() > edge->B.x()){
    return edge->A;
  }
  return edge->B;
}

/*
 * I never sorted the edges... is that bad?
 * I was certain it wouldn't work unless i sorted edges, but it works 
 * I guess there is an inherent sort somewhere in my algorithm
*/
void sortEdges(Edge* edges [], int count){

}



/* returns inverted form of matrix */
float** invertMatrix(float** matrix){

  int rows = 3;
  int columns = 3;

  float** m = createIdentityMatrix(rows, columns);
  
  //hard code the elements, then multiply by determinant in loop
  m[0][0] = matrix[1][1];
  m[0][1] = -matrix[0][1];
  m[0][2] = matrix[0][1]*matrix[1][2] - matrix[1][1]*matrix[0][2];

  m[1][0] = -matrix[1][0];
  m[1][1] = matrix[0][0];
  m[1][2] = matrix[1][0]*matrix[0][2] - matrix[0][0]*matrix[1][2];

  float determinant = matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];

  for(int r = 0; r < rows; r++){
    for(int c = 0; c < columns; c++){
      m[r][c] = m[r][c]*(1/determinant);
    } 
  }

  return m;

}

void printMatrix(float** matrix){

  int rows = 3;
  int columns = 3;  

  cout <<"\n"; 
  for(int r = 0; r < rows; r++){
    for(int c = 0; c < columns; c++){
      cout << matrix[r][c]; 
      cout << " ";
    }
      cout <<"\n"; 
  }

  cout <<"\n"; 
}

/* Multiplies 3x3 matricies A and B 
 * Only for coordinate transform matricies (i.e. bottom row [0,0,1])
*/
float** multiplyMatricies(float** A, float** B){
  float** m = createIdentityMatrix(3,3);
  // because its only 6 elements, we use the analytical formula 
  m[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];  
  m[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1]; 
  m[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]; 


  m[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];  
  m[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1]; 
  m[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]; 

  return m;


}

int clamp(int min, int value, int max){
  return std::max(min, std::min(value, max));
}


float** createScaleMatrix(const GBitmap& bitmap, const GRect rect){

  float** transform = createIdentityMatrix(3, 3);
  
  float Sx = rect.width()/bitmap.width();
  float Sy = rect.height()/bitmap.height();

  transform[0][0] = Sx;
  transform[1][1] = Sy; 

  return transform;

}


float** createSrcDst(const GBitmap& bitmap, const GRect& rect){
 
  float** transform = createIdentityMatrix(3, 3);
  
  float Sx = rect.width()/bitmap.width();
  float Sy = rect.height()/bitmap.height();

  float Tx = 0-rect.left();
  float Ty = 0-rect.top();

  transform[0][0] = Sx;
  transform[1][1] = Sy;

  transform[0][2] = -Tx; 
  transform[1][2] = -Ty; 

  return transform;

}

float getDstX(float x, float y, float** transform){
  float x_coord = x*transform[0][0] + y*transform[0][1] + transform[0][2];
  return x_coord; 
}

float getDstY(float x, float y, float** transform){
  float y_coord = x*transform[1][0] + y*transform[1][1] + transform[1][2];
  return y_coord; 
}

void MyCanvas::fillConvexPolygonBitmap(Edge* edges [], int count, const GBitmap& src, float** untransform){

  //Find highest and lowest points
  GPoint maxYPt = GPoint::Make(the_rect.left(), the_rect.bottom()); // closest to top
  GPoint minYPt = GPoint::Make(the_rect.left(), the_rect.top()); // closest to bottom
 
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(getUpperPoint(edges[i]).y() <= maxYPt.y()){ 
        maxYPt = getUpperPoint(edges[i]);
      }
      if(getLowerPoint(edges[i]).y() >= minYPt.y()){
        minYPt = getLowerPoint(edges[i]);
      }
    }
  }
 
  // these will be the two highest hedges
  Edge* first_edge = new Edge;
  Edge* second_edge = new Edge;

  int highestEdgeIndex = -1;
  // get highest edge
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(getUpperPoint(edges[i]).y() <= maxYPt.y()){  
        first_edge = edges[i];
        highestEdgeIndex = i;
      }
    }
  }

  //second edge = minedge
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(getLowerPoint(edges[i]).y() >= minYPt.y()){  
        second_edge = edges[i];
      }
    }
  }
  
  // get second highest edge
  for(int i = 0; i < count*3; i++){
    if(edges[i] != NULL){
      if(getUpperPoint(edges[i]).y() <= getUpperPoint(second_edge).y()){        
        if(i != highestEdgeIndex){  
          second_edge = edges[i];
        }
      }
    }
  }

  //see which is left edge, which is right edge
  Edge* left_edge = new Edge;
  Edge* right_edge = new Edge;

  if(getLowerPoint(first_edge).x() < getLowerPoint(second_edge).x()){
    left_edge = first_edge;
    right_edge = second_edge;
  }else{
    left_edge = second_edge;
    right_edge = first_edge;    
  }

  for(int y = maxYPt.y(); y < minYPt.y(); y++){
  
    int leftBoundary = computeEdgeBoundary(left_edge, y);
    int rightBoundary = computeEdgeBoundary(right_edge, y);
  
    GPixel* dst = the_bitmap.fPixels;
    dst = (GPixel*)((char*)dst + the_bitmap.rowBytes()*y); 

    for(int x = leftBoundary; x < rightBoundary; x++){ 

      
      GPixel src_pixel = lookupPixel(x,y, untransform, src);
      GPixel dst_pixel = dst[x];
      dst[x] = srcOverBlend(src_pixel, dst_pixel);
    }  
    
    if( y >= getLowerPoint(left_edge).y()){ 
      left_edge = getBoundaryEdge(edges, count, left_edge);
    }
    if( y >= getLowerPoint(right_edge).y()){
      right_edge = getBoundaryEdge(edges, count, right_edge);
    } 
  }
}

class BitmapShader : public GShader{
  public:
  BitmapShader(const GBitmap& bitmap, const float localMatrix[6]){

    float** m = createIdentityMatrix(3,3);
    m[0][0] = localMatrix[0];    
    m[0][1] = localMatrix[1];    
    m[0][2] = localMatrix[2];    
    m[1][0] = localMatrix[3];    
    m[1][1] = localMatrix[4];    
    m[1][2] = localMatrix[5];    

    this->localMatrix = m;
    this->the_bitmap = bitmap;
    this->cm = createIdentityMatrix(3,3);

  }
  bool setContext(const float mat[6]) override {
    
    float** m = createIdentityMatrix(3,3);
    m[0][0] = mat[0];    
    m[0][1] = mat[1];    
    m[0][2] = mat[2];    
    m[1][0] = mat[3];    
    m[1][1] = mat[4];    
    m[1][2] = mat[5];    
    this->cm = m; 

    return true;
  }  
  void shadeRow(int x, int y, int count, GPixel row[]) override { 
    
    for (int i = 0; i < count; ++i) {      
      float** untransform =  invertMatrix(multiplyMatricies(this->cm, this->localMatrix));
      GPixel yourPixel = lookupPixel(x,y, untransform, the_bitmap);         
      row[i] = yourPixel; 
      x++;
    } 
  }
  protected:
    float** cm;
    float** localMatrix;
    GBitmap the_bitmap; 
};
GShader* GShader::FromBitmap(const GBitmap& bitmap, const float localMatrix[6]){
  return new BitmapShader(bitmap, localMatrix);
}
    
class RadialGradientShader : public GShader{
  public:
  RadialGradientShader(const GPoint& center, float radius, const GColor colors[2]){
    this->center = center;
    this->radius = radius;
  
    
   

    this->colors[0] = colors[0];
    this->colors[1] = colors[1];
    this->ctm = createIdentityMatrix(3,3);     
  } 
  bool setContext(const float mat[6]) override { 
    float** m = createIdentityMatrix(3,3);
    m[0][0] = mat[0];    
    m[0][1] = mat[1];    
    m[0][2] = mat[2];    
    m[1][0] = mat[3];    
    m[1][1] = mat[4];    
    m[1][2] = mat[5];    
    this->ctm = m; 
     
    return true;
  }  

  void shadeRow(int x, int y, int count, GPixel row[]) override {
    for (int i = 0; i < count; ++i) {


      int cX =  getDstX(this->center.x(), this->center.y(),this->ctm); 
      int cY =  getDstY(this->center.x(), this->center.y(),this->ctm); 

      GPoint centerCtm = GPoint::Make(cX,cY);
      GPoint kpt = GPoint::Make(x,y);

      int eX = this->center.x();
      int eY = this->center.y() + this->radius;
      GPoint edgePoint = GPoint::Make(eX, eY); 

      int eTX = getDstX(edgePoint.x(), edgePoint.y(), this->ctm);
      int eTY = getDstY(edgePoint.x(), edgePoint.y(), this->ctm);

      GPoint ctmEdgePoint = GPoint::Make(eTX, eTY); 

      float unitValue =  distanceBetween(kpt, centerCtm)/distanceBetween(centerCtm, ctmEdgePoint); 

      if(unitValue > 1){
        unitValue = 1;
      }
           
      GColor combinedColor = linearCombineColors(unitValue, this->colors);
      GPixel yourPixel = makePixelFromColor(combinedColor);
      row[i] = yourPixel; 
      x++;
    }
  }
  protected:
    float** ctm;
    GPoint center;
    float radius;
    GColor colors [2];
};

GShader* GShader::FromRadialGradient(const GPoint& center, float radius, const GColor colors[2]){
  return new RadialGradientShader(center, radius, colors);
}


class LinearGradientShader : public GShader{
  public:
  LinearGradientShader(const GPoint pts[2], const GColor colors[2]){

    this->pts[0] = pts[0];
    this->pts[1] = pts[1];

    this->colors[0] = colors[0];
    this->colors[1] = colors[1];
     

    this->ctm = createIdentityMatrix(3,3);    
  }  
  bool setContext(const float mat[6]) override { 
    float** m = createIdentityMatrix(3,3);
    m[0][0] = mat[0];    
    m[0][1] = mat[1];    
    m[0][2] = mat[2];    
    m[1][0] = mat[3];    
    m[1][1] = mat[4];    
    m[1][2] = mat[5];    
    this->ctm = m; 
     
    return true;
  }  
  void shadeRow(int x, int y, int count, GPixel row[]) override {  
 

    int x0 =  getDstX(pts[0].x(), pts[0].y(),this->ctm); 
    int y0 =  getDstY(pts[0].x(), pts[0].y(),this->ctm); 

    int x1 =  getDstX(pts[1].x(), pts[1].y(),this->ctm); 
    int y1 =  getDstY(pts[1].x(), pts[1].y(),this->ctm); 

    GPoint ctmPts[2];
    ctmPts[0] = GPoint::Make(x0,y0);
    ctmPts[1] = GPoint::Make(x1,y1);

    float** M = mapToLine(ctmPts);
    for (int i = 0; i < count; ++i) {
      GPoint kpt = GPoint::Make(x,y);
      GPoint mappedPt = projectPointontoLine(ctmPts, kpt);                     
      float unitValue = getUnitValue(mappedPt.x(), mappedPt.y(), M);     
      GColor combinedColor = linearCombineColors(unitValue, this->colors);
      GPixel yourPixel = makePixelFromColor(combinedColor);
      row[i] = yourPixel;      
      x++;
    }
  }

  protected:
    float** ctm;
    GPoint pts [2];
    GColor colors [2]; 
};

GShader* GShader::FromLinearGradient(const GPoint pts[2], const GColor colors[2]){
  return new LinearGradientShader( pts, colors); // a LinearGradientShader object
}

float** mapToLine(const GPoint* pts){

  float** M = createIdentityMatrix(3,3); 

  float dx = pts[1].x() - pts[0].x();
  float dy = pts[1].y() - pts[0].y();

  M[0][0] = dx;
  M[0][1] = -dy;
  M[0][2] = pts[0].x();
  M[1][0] = dy; 
  M[1][1] = dx;
  M[1][2] = pts[0].y();


  return M;

}

float getUnitValue(int x, int y, float** M){
  float** transform = invertMatrix(M);
  float x_coord = x*transform[0][0] + y*transform[0][1] + transform[0][2];
  float y_coord = x*transform[1][0] + y*transform[1][1] + transform[1][2];
  float z_coord = x*transform[2][0] + y*transform[2][1] + transform[2][2];
  return x_coord; 
}

GPoint projectPointontoLine(const GPoint* pts, GPoint kpt){

  GPoint p0 = pts[0];
  GPoint p1 = pts[1];

  float adx = kpt.x() - p0.x();
  float ady = kpt.y() - p0.y();
  float alength = sqrt( adx*adx + ady*ady);
  if(alength == 0){ return pts[0]; }

  float bdx = p1.x() - p0.x();
  float bdy = p1.y() - p0.y();
  float blength = sqrt( bdx*bdx + bdy*bdy);
 
  float adotb = adx*bdx + ady*bdy;
  float costheta = (adotb)/(alength*blength);
 
  float projectedx = (bdx/blength)*costheta*alength; 
  float projectedy = (bdy/blength)*costheta*alength;
  
  float projectedlength = sqrt( projectedx*projectedx + projectedy*projectedy);

  if(costheta < 0){
    return pts[0]; 
  }else if(projectedlength > blength){ 
    return pts[1];
  }
  return GPoint::Make(projectedx + p0.x(), projectedy + p0.y());
}

GColor linearCombineColors(float unitValue, const GColor* colors){

  // linear combine 2 colors based on unit value -> colorA*(1-unitvalue) + colorB*unitvalue 
  // if unitvalue is 1, color is all color B
  GColor pinnedColors[2];

  pinnedColors[0] = colors[0];//.pinToUnit();
  pinnedColors[1] = colors[1];//.pinToUnit();

  if(unitValue == 0) {return pinnedColors[0];} 
  if(unitValue == 1) {return pinnedColors[1];}

  float lcA = pinnedColors[0].fA*(1-unitValue) + pinnedColors[1].fA*(unitValue);
  float lcR = pinnedColors[0].fR*(1-unitValue) + pinnedColors[1].fR*(unitValue);
  float lcG = pinnedColors[0].fG*(1-unitValue) + pinnedColors[1].fG*(unitValue);
  float lcB = pinnedColors[0].fB*(1-unitValue) + pinnedColors[1].fB*(unitValue);
 
 
  return GColor::MakeARGB(lcA, lcR, lcG, lcB);
}

// returns distance between A and B 
float distanceBetween(GPoint A, GPoint B){
  
  float dx = A.x() - B.x();
  float dy = A.y() - B.y();
  float dst = sqrt( dx*dx + dy*dy);
  return dst;
}


/**
 *  Stroke the specified polygon using the Stroke settings. If isClosed is true, then the
 *  drawn stroke should connect the first and last points of the polygon, else it should not,
 *  and those end-caps should reflect the Stroke.fAddCap setting.
 */
void MyCanvas::strokePolygon(const GPoint pts [], int count, bool isClosed, const Stroke& stroke, GShader* shader){

  if ((isClosed && count < 3) || count < 2)
    return; 




  Edge** edges = new Edge*[count];
  
  //fillConvexPolygon(pts, count, GColor::MakeARGB(1, 0, 0, 0));

  for(int i = 0; i < count; i++){
    edges[i] = NULL;
  } 

  for(int i = 0; i < count; i++){
    Edge* edge = new Edge;
    edge->A = pts[i];
    edge->B = pts[(i+1) % count]; 
    edges[i] = edge;
  } 



  GPoint qpts[4]; // the points of the rect we made from edge
  GPoint bpts[4]; // for getting the pts for the bevel (never drawn)
  GPoint tpts[3]; // for bevel triangle notch
  GPoint cpts[4]; // for caps
  GPoint mpts[4]; // for meiters

  for(int i = 0; i < count; i++){ // there's only count edges, we didn't clipp anything

    Edge* edge = edges[i%count]; 
    GPoint A = edge->A;
    GPoint B = edge->B;

    float x = B.x() - A.x(); 
    float y = B.y() - A.y();
    
    GPoint u = GPoint::Make(-x,-y);

    u = makeUnitVector(u);


    float len = sqrt(x*x + y*y);
    float rad = stroke.fWidth/2 ;

    float xprime = x*rad/len;
    float yprime = y*rad/len;
  
    GPoint abt = GPoint::Make(-yprime, xprime);
   
    GPoint alpha = GPoint::Make(abt.x(), abt.y());
    
    qpts[0] = GPoint::Make(A.x() + abt.x(), A.y() + abt.y());
    qpts[1] = GPoint::Make(A.x() - abt.x(), A.y() - abt.y());
    qpts[2] = GPoint::Make(B.x() - abt.x(), B.y() - abt.y());
    qpts[3] = GPoint::Make(B.x() + abt.x(), B.y() + abt.y()); 
    
    
     
    
    if(i == 0 || i == count-1){  
    if(stroke.fAddCap){
    if(!isClosed){ 
        cpts[0] = GPoint::Make(qpts[0].x()-xprime, qpts[0].y()-yprime);   
        cpts[1] = GPoint::Make(qpts[1].x()-xprime, qpts[1].y()-yprime);
        cpts[2] = qpts[1];
        cpts[3] = qpts[0];   
        //fillConvexPolygon(cpts, 4, GColor::MakeARGB(1, 0, 0, 0));
        shadeConvexPolygon(cpts, 4, shader);
    }
    }
    } 
    
    



    if(!isClosed && i == count-1){return;}
    
    
    // get point from next edge    
    edge = edges[(i+1)%count]; 
    A = edge->A;
    B = edge->B;
    
    x = B.x() - A.x(); 
    y = B.y() - A.y();

    GPoint v = GPoint::Make(x,y);
    v = makeUnitVector(v);

    len = sqrt(x*x + y*y);
    rad = stroke.fWidth /2;

    xprime = x*rad/len;
    yprime = y*rad/len;
  
    abt = GPoint::Make(-yprime, xprime);
    
    GPoint beta = GPoint::Make(abt.x(), abt.y());

    bpts[0] = GPoint::Make(A.x() + abt.x(), A.y() + abt.y());
    bpts[1] = GPoint::Make(A.x() - abt.x(), A.y() - abt.y());
    bpts[2] = GPoint::Make(B.x() - abt.x(), B.y() - abt.y());
    bpts[3] = GPoint::Make(B.x() + abt.x(), B.y() + abt.y()); 
 
    tpts[0] = A;
    tpts[2] = qpts[2]; 
    tpts[1] = bpts[1];



   
    shadeConvexPolygon(qpts, 4, shader);
    //fillConvexPolygon(tpts, 3, GColor::MakeARGB(1, 0, 0, 0));


    GPoint vA = vectorFromEdge(edges[i%count]);
    GPoint vB = vectorFromEdge(edges[(i+1)%count]);

    //float acrossb = vA.x()*vB.y() - vA.y()*vB.x();
    float acrossb = u.x()*v.y() - u.y()*v.x();
    if(acrossb>0){//invert u
      alpha = GPoint::Make(-alpha.x(), -alpha.y());
      beta = GPoint::Make(-beta.x(), -beta.y()); 
    }else{//invert v

    } 

    float ml = getMiterLength(u, v, stroke.fWidth, stroke.fMiterLimit);


    if(ml == -1){
      shadeConvexPolygon(tpts, 3, shader); 
    }else{
      // make some meiter points and shade them;

        

      

      GPoint mdirection = GPoint::Make(-alpha.x()-beta.x(),-alpha.y()-beta.y());
      GPoint munit = makeUnitVector(mdirection);
      GPoint miterTip = GPoint::Make(munit.x()*ml, munit.y()*ml); 
      mpts[1] = A;
      mpts[2] = qpts[2]; 
      mpts[3] = GPoint::Make(miterTip.x()+A.x(), miterTip.y()+A.y());
      mpts[0] = bpts[1];
 
      //fillConvexPolygon(mpts, 4, GColor::MakeARGB(1, 0, 0, 0));
      shadeConvexPolygon(mpts, 4, shader);
    }  



  }  
  
  return;
} 

//Returns length, as if A is vector from origin
float getVectorLength(GPoint A){
 
  GPoint B = GPoint::Make(0,0);
 
  float dx = A.x() - B.x();
  float dy = A.y() - B.y(); 
  float dst = sqrt( dx*dx + dy*dy);

  return dst;
}

GPoint translatePoint(GPoint pt, float** matrix){
  GPoint tpt = GPoint::Make(getDstX(pt.x(), pt.y(), matrix), getDstY(pt.x(), pt.y(), matrix)); 
  return tpt;
}

GPoint makeUnitVector(GPoint v){

  return GPoint::Make(v.x()/getVectorLength(v), v.y()/getVectorLength(v));

}

float getMiterLength(GPoint u, GPoint v, float width, float limit){

  float h = sqrt(2/ (1 - ((u.x()*v.x()+u.y()*v.y()))));
  if(h > limit){
    return -1;
  }

  return (width/2) * h;

}


GPoint vectorFromEdge(Edge* e){

  float x = e->B.x() - e->A.x();
  float y = e->B.y() - e->A.y();

  return GPoint::Make(x,y);

}









