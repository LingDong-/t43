#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define N_PIXEL_NEIGHBOR 8

typedef struct _point_t {
  float x;
  float y;
  struct _point_t * next;
  float curv;
  float nx;
  float ny;
} point_t;

typedef struct _polyline_t {
  point_t* head;
  point_t* tail;
  struct _polyline_t* prev;
  struct _polyline_t* next;
  int size;
  int id;
  int parent;
  char is_hole;
} polyline_t;

polyline_t* new_polyline(){
  polyline_t* q0 = (polyline_t*)calloc(1,sizeof(polyline_t));
  return q0;
}

void add_point_to_polyline(polyline_t* q, float x, float y){
  point_t* p = (point_t*)malloc(sizeof(point_t));
  p->x = x;
  p->y = y;
  p->next = NULL;
  if (!q->head){
    q->head = p;
    q->tail = p;
  }else{
    q->tail->next = p;
    q->tail = p;
  }
  q->size++;
}

polyline_t* prepend_polyline(polyline_t* q0, polyline_t* q1){
  if (!q0){
    return q1;
  }
  q1->next = q0;
  q0->prev = q1;
  return q1;
}

void pop_head_polyline(polyline_t* q){
  point_t* p = q->head;
  q->head = q->head->next;
  free(p);
  q->size --;
}

void cat_tail_polyline(polyline_t* q0, polyline_t* q1){
  if (!q1){
    return;
  }
  if (!q0){
    *q0 = *new_polyline();
  }
  if (!q0->head){
    q0->head = q1->head;
    q0->tail = q1->tail;
    return;
  }
  q0->tail->next = q1->head;
  q0->tail = q1->tail;
  q0->size += q1->size;
  q0->tail->next = NULL;
}

void reverse_polyline(polyline_t* q){
  if (!q || (q->size < 2)){
    return;
  }
  if (q->size < 3){
    point_t* q_head = q->head;
    q->tail->next = q_head;
    q->head = q->tail;
    q->tail = q_head;
    q->tail->next = NULL;
    return;
  }
  q->tail->next = q->head;
  point_t* it0 = q->head;
  point_t* it1 = it0->next;
  point_t* it2 = it1->next;
  int i; for (i = 0; i < q->size-1; i++){
    it1->next = it0;
    it0 = it1;
    it1 = it2;
    it2 = it2->next;
  }
  point_t* q_head = q->head;
  q->head = q->tail;
  q->tail = q_head;
  q->tail->next = NULL;
}

polyline_t* dup_polyline(polyline_t* p){
  polyline_t*q = new_polyline();
  point_t* it = p->head;
  while (it){
    add_point_to_polyline(q,it->x,it->y);
    it = it->next;
  }
  return q;
}

void destroy_polylines(polyline_t* q){
  if (!q){
    return;
  }
  polyline_t* it = q;
  while(it){
    polyline_t* lt = it->next;
    point_t* jt = it->head;
    while(jt){
      point_t* kt = jt->next;
      free(jt);
      jt = kt;
    }
    free(it);
    it = lt;
  }
}
void polylines_to_svg(const char* pth, int xmin, int ymin, int xmax, int ymax, polyline_t* q){
  FILE *fp;
  fp = fopen(pth, "w+");
  int w = xmax-xmin;
  int h = ymax-ymin;
  fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"%d %d %d %d\">\n",w,h,xmin,ymin,w,h);

  polyline_t* it = q;
  while(it){
    fprintf(fp, "<path fill=\"none\" stroke-width=\"1\" stroke=\"#%06x\" stroke-linecap=\"round\" stroke-linejoin=\"round\" d=\"M",rand()%0xFFFFFF);
    point_t* jt = it->head;
    while(jt){
      fprintf(fp,"%f,%f ",jt->x,jt->y);
      jt = jt->next;
    }
    fprintf(fp, "\" />\n");
    it = it->next;
    // break;
  }
  fprintf(fp, "</svg>\n");
  fclose(fp);
}

int neighbor_id_to_index(int i, int j, int id, int* oi, int* oj){
  if (id == 0){return *oi = i,   *oj = j+1;}
  if (id == 1){return *oi = i-1, *oj = j+1;}
  if (id == 2){return *oi = i-1, *oj = j;  }
  if (id == 3){return *oi = i-1, *oj = j-1;}
  if (id == 4){return *oi = i,   *oj = j-1;}
  if (id == 5){return *oi = i+1, *oj = j-1;}
  if (id == 6){return *oi = i+1, *oj = j;  }
  if (id == 7){return *oi = i+1, *oj = j+1;}
  return *oi=-1, *oj=-1;
}
int neighbor_index_to_id(int i0, int j0, int i, int j){
  int di = i - i0;
  int dj = j - j0;
  if (di == 0 && dj == 1){return 0;}
  if (di ==-1 && dj == 1){return 1;}
  if (di ==-1 && dj == 0){return 2;}
  if (di ==-1 && dj ==-1){return 3;}
  if (di == 0 && dj ==-1){return 4;}
  if (di == 1 && dj ==-1){return 5;}
  if (di == 1 && dj == 0){return 6;}
  if (di == 1 && dj == 1){return 7;}
  return -1;
}

int ccw_non0(int8_t* F, int w, int h, int i0, int j0, int i, int j, int offset, int *oi, int *oj){
  int id = neighbor_index_to_id(i0,j0,i,j);
  for (int k = 0; k < N_PIXEL_NEIGHBOR; k++){
    int kk = (k+id+offset + N_PIXEL_NEIGHBOR*2) % N_PIXEL_NEIGHBOR;
    int i, j;
    neighbor_id_to_index(i0,j0,kk,&i,&j);
    if (F[i*w+j]!=0){
      return *oi = i, *oj = j;
    }
  }
  return *oi=-1, *oj=-1;
}

int cw_non0(int8_t* F, int w, int h, int i0, int j0, int i, int j, int offset, int *oi, int *oj){
  int id = neighbor_index_to_id(i0,j0,i,j);
  for (int k = 0; k < N_PIXEL_NEIGHBOR; k++){
    int kk = (-k+id-offset + N_PIXEL_NEIGHBOR*2) % N_PIXEL_NEIGHBOR;
    int i, j;
    neighbor_id_to_index(i0,j0,kk,&i,&j);
    if (F[i*w+j]!=0){
      return *oi = i, *oj = j;
    }
  }
  return *oi=-1, *oj=-1;
}

polyline_t* find_contours(int8_t* F, int w, int h) {
  // Topological Structural Analysis of Digitized Binary Images by Border Following.
  // Suzuki, S. and Abe, K., CVGIP 30 1, pp 32-46 (1985)
  int nbd = 1;
  int lnbd = 1;
  
  polyline_t* contours = NULL;
  
  // Without loss of generality, we assume that 0-pixels fill the frame
  // of a binary picture
  for (int i = 1; i < h-1; i++){
    F[i*w] = 0; F[i*w+w-1]=0;
  }
  for (int i = 0; i < w; i++){
    F[i] = 0; F[w*h-1-i]=0;
  }
  
  //Scan the picture with a TV raster and perform the following steps
  //for each pixel such that fij # 0. Every time we begin to scan a
  //new row of the picture, reset LNBD to 1.
  for (int i = 1; i < h-1; i++) {

    lnbd = 1;
    
    for (int j = 1; j < w-1; j++) {
      
      int i2 = 0, j2 = 0;
      if (F[i*w+j] == 0) {
        continue;
      }
      //(a) If fij = 1 and fi, j-1 = 0, then decide that the pixel
      //(i, j) is the border following starting point of an outer
      //border, increment NBD, and (i2, j2) <- (i, j - 1).
      if (F[i*w+j] == 1 && F[i*w+(j-1)] == 0) {
        nbd ++;
        i2 = i;
        j2 = j-1;
        
        
        //(b) Else if fij >= 1 and fi,j+1 = 0, then decide that the
        //pixel (i, j) is the border following starting point of a
        //hole border, increment NBD, (i2, j2) <- (i, j + 1), and
        //LNBD + fij in case fij > 1.
      } else if (F[i*w+j]>=1 && F[i*w+j+1] == 0) {
        nbd ++;
        i2 = i;
        j2 = j+1;
        if (F[i*w+j]>1) {
          lnbd = F[i*w+j];
        }
        
        
      } else {
        //(c) Otherwise, go to (4).
        //(4) If fij != 1, then LNBD <- |fij| and resume the raster
        //scan from pixel (i,j+1). The algorithm terminates when the
        //scan reaches the lower right corner of the picture
        if (F[i*w+j]!=1){lnbd = abs(F[i*w+j]);}
        continue;
        
      }
      //(2) Depending on the types of the newly found border
      //and the border with the sequential number LNBD
      //(i.e., the last border met on the current row),
      //decide the parent of the current border as shown in Table 1.
      // TABLE 1
      // Decision Rule for the Parent Border of the Newly Found Border B
      // ----------------------------------------------------------------
      // Type of border B'
      // \    with the sequential
      //     \     number LNBD
      // Type of B \                Outer border         Hole border
      // ---------------------------------------------------------------
      // Outer border               The parent border    The border B'
      //                            of the border B'
      //
      // Hole border                The border B'      The parent border
      //                                               of the border B'
      // ----------------------------------------------------------------
      
      polyline_t* B = new_polyline();

      add_point_to_polyline(B,j,i);
      
      B->is_hole = (j2 == j+1);
      B->id = nbd;
      contours = prepend_polyline(contours,B);
      
      
      polyline_t* B_0 = NULL;
      polyline_t* it = contours;
      while (it){
        if (it->id == lnbd){
          B_0 = it;
          break;
        }
        it = it->next;
      }
      if (B_0){
        if (B_0->is_hole){
          if (B->is_hole){
            B->parent = B_0->parent;
          }else{
            B->parent = lnbd;
          }
        }else{
          if (B->is_hole){
            B->parent = lnbd;
          }else{
            B->parent = B_0->parent;
          }
        }
      }else{
        if (B->is_hole){
          B->parent = lnbd;
        }else{
          B->parent = 0;
        }
      }
      
      //(3) From the starting point (i, j), follow the detected border:
      //this is done by the following substeps (3.1) through (3.5).
      
      //(3.1) Starting from (i2, j2), look around clockwise the pixels
      //in the neigh- borhood of (i, j) and tind a nonzero pixel.
      //Let (i1, j1) be the first found nonzero pixel. If no nonzero
      //pixel is found, assign -NBD to fij and go to (4).
      int i1 = -1, j1 = -1;
      cw_non0(F,w,h,i,j,i2,j2,0,&i1,&j1);
      if (i1 == -1){
        F[i*w+j] = -nbd;
        //go to (4)
        if (F[i*w+j]!=1){lnbd = abs(F[i*w+j]);}
        continue;
      }

      // (3.2) (i2, j2) <- (i1, j1) ad (i3,j3) <- (i, j).
      i2 = i1;
      j2 = j1;
      int i3 = i;
      int j3 = j;
      
      
      while (1){
        //(3.3) Starting from the next elementof the pixel (i2, j2)
        //in the counterclock- wise order, examine counterclockwise
        //the pixels in the neighborhood of the current pixel (i3, j3)
        //to find a nonzero pixel and let the first one be (i4, j4).
        
        int i4, j4;
        ccw_non0(F,w,h,i3,j3,i2,j2,1,&i4,&j4);
        
        add_point_to_polyline(contours, j4,i4);
        
        //(a) If the pixel (i3, j3 + 1) is a O-pixel examined in the
        //substep (3.3) then fi3, j3 <-  -NBD.
        if (F[i3*w+j3+1] == 0){
          F[i3*w+j3] = -nbd;
          
          //(b) If the pixel (i3, j3 + 1) is not a O-pixel examined
          //in the substep (3.3) and fi3,j3 = 1, then fi3,j3 <- NBD.
        }else if (F[i3*w+j3] == 1){
          F[i3*w+j3] = nbd;
        }else{
          //(c) Otherwise, do not change fi3, j3.
        }
        
        //(3.5) If (i4, j4) = (i, j) and (i3, j3) = (i1, j1)
        //(coming back to the starting point), then go to (4);
        if (i4 == i && j4 == j && i3 == i1 && j3 == j1){
          if (F[i*w+j]!=1){lnbd = abs(F[i*w+j]);}
          break;
          
          //otherwise, (i2, j2) + (i3, j3),(i3, j3) + (i4, j4),
          //and go back to (3.3).
        }else{
          i2 = i3;
          j2 = j3;
          i3 = i4;
          j3 = j4;
        }
      }
    }
  }
  return contours;
}


float point_distance_to_segment(point_t p, point_t p0, point_t p1) {
  // https://stackoverflow.com/a/6853926
  float x = p.x;   float y = p.y;
  float x1 = p0.x; float y1 = p0.y;
  float x2 = p1.x; float y2 = p1.y;
  float A = x - x1; float B = y - y1; float C = x2 - x1; float D = y2 - y1;
  float dot = A*C+B*D;
  float len_sq = C*C+D*D;
  float param = -1;
  if (len_sq != 0) {
    param = dot / len_sq;
  }
  float xx; float yy;
  if (param < 0) {
    xx = x1; yy = y1;
  }else if (param > 1) {
    xx = x2; yy = y2;
  }else {
    xx = x1 + param*C;
    yy = y1 + param*D;
  }
  float dx = x - xx;
  float dy = y - yy;
  return sqrt(dx*dx+dy*dy);
}

void approx_poly_dp(polyline_t* polyline, float epsilon){
  // https://en.wikipedia.org/wiki/Ramer–Douglas–Peucker_algorithm

  if (polyline->size < 3){
    return;
  }
  float dmax   = 0;
  point_t* argmax;
  point_t* argmax_prev;
  int imax = 0;
  int i = 0;
  point_t* it_prev = polyline->head;
  point_t* it = polyline->head->next;
  while (it->next){
    float d = point_distance_to_segment(*it,*(polyline->head),*(polyline->tail));
    if (d > dmax){
      dmax = d;
      argmax = it;
      argmax_prev = it_prev;
      imax = i;
    }
    it = it->next;
    it_prev = it_prev->next;
    i++;
  }
  if (dmax > epsilon){
    polyline_t R;
    R.head = argmax;
    R.tail = polyline->tail;
    R.size = polyline->size-imax;

    polyline->tail = argmax_prev;
    argmax_prev->next = NULL;
    polyline->size = imax;

    add_point_to_polyline(polyline, argmax->x, argmax->y);

    approx_poly_dp(polyline,epsilon);
    approx_poly_dp(&R,epsilon);
    pop_head_polyline(&R);
    cat_tail_polyline(polyline,&R);
    
  }else{
    it = polyline->head->next;
    while (it->next){
      point_t* nxt = it->next;
      free(it);
      it = nxt;
    }
    polyline->head->next = polyline->tail;
    polyline->size = 2;
  }
}


int isect_circ_line(float cx,float cy,float r,float x0,float y0,float x1,float y1,float*tout){
  //https://stackoverflow.com/a/1084899
  float dx = x1-x0;
  float dy = y1-y0;
  float fx = x0-cx;
  float fy = y0-cy;
  float a = dx*dx+dy*dy;
  float b = 2*(fx*dx+fy*dy);
  float c = (fx*fx+fy*fy)-r*r;
  float discriminant = b*b-4*a*c;
  if (discriminant<0){
    return 0;
  }
  discriminant = sqrt(discriminant);
  float t0 = (-b - discriminant)/(2*a);
  if (0 <= t0 && t0 <= 1){
    return 0;
  }
  float t = (-b + discriminant)/(2*a);
  if (t > 1 || t < 0){
    return 0;
  }
  (*tout) = t;
  return 1;
}

void resample(polyline_t* polyline,float step){
  if (polyline->size < 2){
    return;
  }

  polyline_t* out = new_polyline();
  add_point_to_polyline(out, polyline->head->x, polyline->head->y);

  point_t* next;
  point_t* it = polyline->head;
  while(it->next){
    point_t* a = it;
    point_t* b = it->next;
    float dx = b->x-a->x;
    float dy = b->y-a->y;
    float d = sqrt(dx*dx+dy*dy);
    if (d == 0){
      it = it->next;
      continue;
    }
    int n = (d/step);
    float rest = (n*step)/d;
    float rpx = a->x * (1-rest) + b->x * rest;
    float rpy = a->y * (1-rest) + b->y * rest;
    for (int j = 1; j <= n; j++){
      float t = (float)j/(float)n;
      float x = a->x*(1-t) + rpx*t;
      float y = a->y*(1-t) + rpy*t;
      add_point_to_polyline(out,x,y);
    }
    next = NULL;
    point_t* jt = it->next;
    while(jt->next){
      point_t* b = jt;
      point_t* c = jt->next;
      if (b->x == c->x && b->y == c->y){
        jt = jt->next;
        continue;
      }
      float t;
      int hit = isect_circ_line(rpx,rpy,step,b->x,b->y,c->x,c->y,&t);
      if (!hit){
        jt = jt->next;
        continue;
      }
      
      float qx = b->x*(1-t)+c->x*t;
      float qy = b->y*(1-t)+c->y*t;
      add_point_to_polyline(out,qx,qy);
      jt->x = qx;
      jt->y = qy;
      next = jt;
      break;
    }
    if (next == NULL){
      break;
    }
    it = next;
  }

  float mx = polyline->tail->x;
  float my = polyline->tail->y;
  if (out->size > 1){
    float lx = out->tail->x;
    float ly = out->tail->y;
    float d = sqrt((mx-lx)*(mx-lx)+(my-ly)*(my-ly));
    if (d < step*0.5){
      out->tail->x = mx;
      out->tail->y = my;
    }else{
      add_point_to_polyline(out, mx, my);
    }
  }else{
    add_point_to_polyline(out, mx, my);
  }
  it = polyline->head;
  while (it){
    point_t* next = it->next;
    free(it);
    it = next;
  }
  polyline->size = out->size;
  polyline->head = out->head;
  polyline->tail = out->tail;
  free(out);
}




typedef struct {
  uint8_t* data;
  int w;
  int h;
  int d;
  uint64_t cap;
} rast_t;

rast_t new_rast(int w, int h, int d){
  rast_t r;
  r.cap = (uint64_t)w*(uint64_t)h*(uint64_t)d;
  r.data = (uint8_t*) calloc(r.cap,sizeof(uint8_t));
  r.w = w;
  r.h = h;
  r.d = d;
  return r;
}
inline uint8_t rast_get(rast_t* r, int x, int y, int z){
  if (x < 0 || x >= r->w || y < 0 || y >= r->h || z < 0 || z >= r->d){
    return 0;
  }
  return r->data[z*(r->w*r->h)+y*(r->w)+x];
}
uint8_t rast_get_cp_border(rast_t* r, int x, int y, int z){
  x = MIN(MAX(x,0),r->w-1);
  y = MIN(MAX(y,0),r->h-1);
  z = MIN(MAX(z,0),r->d-1);
  return r->data[z*(r->w*r->h)+y*(r->w)+x];
}
inline uint8_t rast_set(rast_t* r, int x, int y, int z, uint8_t v){
  if (x < 0 || x >= r->w || y < 0 || y >= r->h || z < 0 || z >= r->d){
    return 0;
  }
  return r->data[z*(r->w*r->h)+y*(r->w)+x] = v;
}
void rast_clear(rast_t* r){
  memset(r->data, 0, r->w*r->h*r->d);
}
void rast_destroy(rast_t* r){
  free(r->data);
}
rast_t read_voxels_from_file(const char* pth){
  FILE * fp;
  fp = fopen(pth,"r");
  if (!fp){ printf("cannot open file."); exit(1); }
  int w, h, d;
  fread(&w, 4, 1, fp);
  fread(&h, 4, 1, fp);
  fread(&d, 4, 1, fp);
  int n = w*h*d;
  rast_t r = new_rast(w,h,d);
  fread(r.data,1,n,fp);
  return r;
}

rast_t read_voxels_from_file_with_padding(const char* pth, int padx, int pady){
  FILE * fp;
  fp = fopen(pth,"r");
  if (!fp){ printf("cannot open file."); exit(1); }
  int w, h, d;
  fread(&w, 4, 1, fp);
  fread(&h, 4, 1, fp);
  fread(&d, 4, 1, fp);
  rast_t r = new_rast(w+padx*2,h+pady*2,d);
  for (int i = 0; i < d; i++){
    for (int j = 0; j < h; j++){
      fread(r.data + ( i * (r.w*r.h) + (pady+j)*r.w + padx),1, w, fp);
    }
  }
  return r;
}


void write_voxels_to_file(const char* pth, rast_t* vox){
  FILE *fp;
  if (pth){
    fp = fopen(pth, "wb");
  }else{
    fp = stdout;
  }
  fwrite(&(vox->w),4,1,fp);
  fwrite(&(vox->h),4,1,fp);
  fwrite(&(vox->d),4,1,fp);
  int n = vox->w*vox->h*vox->d;
  fwrite((vox->data),1,n,fp);
  fclose(fp);
}


void write_gif(const char* pth,rast_t* im){
  FILE *fp;
  fp = fopen(pth, "wb");
  fputc(0x47,fp);fputc(0x49,fp);fputc(0x46,fp);fputc(0x38,fp);fputc(0x39,fp);fputc(0x61,fp);
  fputc(im->w&0xff,fp);
  fputc((im->w>>8)&0xff,fp);
  fputc(im->h&0xff,fp);
  fputc((im->h>>8)&0xff,fp);
  fputc(0xf6,fp);
  fputc(0,fp); fputc(0,fp);
  for (int i = 0; i < 127; i++){
    // fputc(i*2,fp);fputc(i*2,fp);fputc(i*2,fp);
    int b = ((i >> 5) & 0b11) * 85;
    int g = ((i >> 2) & 0b111) * 36;
    int r = (i & 0b11) * 85;
    fputc(r,fp);fputc(g,fp);fputc(b,fp);
  }
  fputc(0xff,fp);fputc(0xff,fp);fputc(0xff,fp);
  fputc(0x2c,fp);fputc(0,fp);fputc(0,fp);fputc(0,fp);fputc(0,fp);
  fputc(im->w&0xff,fp);
  fputc((im->w>>8)&0xff,fp);
  fputc(im->h&0xff,fp);
  fputc((im->h>>8)&0xff,fp);
  fputc(0,fp);fputc(7,fp);
  int n = im->w*im->h/126;
  int inc = n*126;
  int exc = im->w*im->h-inc;
  for (int i = 0; i < n; i++){
    fputc(0x7f,fp);
    fputc(0x80,fp);
    for (int j = 0; j < 126; j++){
      // fputc(MIN(127,im->data[i*126+j]),fp);
      fputc(im->data[i*126+j]%128,fp);
    }
  }
  if (exc){
    fputc(exc+1,fp);
    fputc(0x80,fp);
    for (int i = 0; i < exc; i++){
      // fputc(MIN(127,im->data[inc+i]),fp);
      fputc(im->data[inc+i]%128,fp);
    }
  }
  fputc(0x01,fp);fputc(0x81,fp);fputc(0x00,fp);fputc(0x3B,fp);
  fclose(fp);
}

void compute_normal(polyline_t* poly){
  point_t* prev = poly->tail;
  point_t* it = poly->head->next;
  point_t* next = poly->head->next->next;
  while (1){
    float x0 = prev->x;
    float y0 = prev->y;
    float x1 = it->x;
    float y1 = it->y;
    float x2 = next->x;
    float y2 = next->y;

    float a0 = atan2(y0-y1,x0-x1);
    float a1 = atan2(y2-y1,x2-x1);
    float a = (a0+a1)/2;
    while (a < a0){
      a += M_PI*2;
    }
    if (a - M_PI > a0 && a - M_PI - a0 < a - a0){
      a -= M_PI;
    }
    a-=M_PI;
    it->nx = cos(a);
    it->ny = sin(a);
    prev = it;
    it = next;
    next = next->next ? next->next : poly->head->next;
    if (it == poly->head->next){
      break;
    }
  }
  poly->head->curv = poly->tail->curv;
  poly->head->nx = poly->tail->nx;
  poly->head->ny = poly->tail->ny;
}

void compute_curvature(polyline_t* poly){
  point_t* prev = poly->tail;
  point_t* it = poly->head->next;
  point_t* next = poly->head->next->next;
  while (1){
    float x0 = prev->x;
    float y0 = prev->y;
    float x1 = it->x;
    float y1 = it->y;
    float x2 = next->x;
    float y2 = next->y;

    float dx = x1-x0;
    float dy = y1-y0;
    float ex = x2-x1;
    float ey = y2-y1;

    float crs  = dx*ey-dy*ex;
    crs /= sqrt(dx*dx+dy*dy)*sqrt(ex*ex+ey*ey);
    it->curv = -crs;

    prev = it;
    it = next;
    next = next->next ? next->next : poly->head->next;
    if (it == poly->head->next){
      break;
    }
  }
  poly->head->curv = poly->tail->curv;
  poly->head->nx = poly->tail->nx;
  poly->head->ny = poly->tail->ny;
}

void smoothen_curvature(polyline_t* poly,int rad){

  int ksize = rad*2+1;
  float kern[ksize];
  float sigma = 0.3*((float)(ksize-1)*0.5 - 1) + 0.8;
  float ss2 = sigma*sigma*2;
  for (int i = 0; i < ksize; i++){
    float x = (float)(i-ksize/2);
    float z = exp(-(x*x)/(ss2))/(2.5066282746*sigma);
    kern[i] = z;
  }

  float oldcurv[poly->size];
  point_t* it = poly->head;
  int idx = 0;
  while (it){
    oldcurv[idx] = it->curv;
    it = it->next;
    idx++;
  }
  it = poly->head;
  idx = 0;
  while (it){
    float s = 0;
    for (int k = 0; k < ksize; k++){
      int i = idx-rad+k;
      i = (i + poly->size) % poly->size;
      s += kern[k] * oldcurv[i];
    }
    it->curv = s;
    idx ++;
    it = it->next;
  }
}

float poly_area(polyline_t* poly){
  if (poly->size <= 3){
    return 0;
  }
  float a = 0;
  point_t* it = poly->head;
  point_t* next = it->next;
  while (next){
    a += it->x * next->y - next->x * it->y;
    it = next;
    next = next->next;
  }
  return a * 0.5;
}





int seg_isect_y_line(float y, point_t p0, point_t p1, point_t* out){
  if (p0.y > p1.y){
    point_t tmp = p0;
    p0 = p1;
    p1 = tmp;
  }
  if (fabs(p0.y-p1.y)<0.0001){
    p1.y += 0.0001;
  }
  float t = (y - p0.y)/(p1.y-p0.y);
  // if (isnan(t)){
  //   printf("t %f %f %f %f\n",y,p0.y,p1.y,p0.y);
  // }
  if (t < 0 || t >= 1){
    return 0;
  }
  out->x = p0.x * (1-t) + p1.x * t;
  out->y = y;
  out->next = NULL;
  return 1;
}

uint8_t* tmp_data = NULL;
int tmp_data_len = 0;
int tmp_data_cap = 0;

void tmp_data_clear(){
  tmp_data_len = 0;
}

void* tmp_data_add(int size){
  if (!tmp_data){
    tmp_data_cap = 1024;
    tmp_data = (uint8_t*)malloc(1024*sizeof(uint8_t));
    tmp_data_len = 0;
  }
  if (tmp_data_len>=tmp_data_cap){
    tmp_data_cap += 1024;
    tmp_data = (uint8_t*)realloc(tmp_data,tmp_data_cap*sizeof(uint8_t));
  }
  tmp_data_len += size;
  return &tmp_data[tmp_data_len-size];
}

int cmp_ptx (const void * a, const void * b) {
  float x = *(float*)a - *(float*)b ;
  return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}


void draw_polygon_fills(rast_t* im, polyline_t* polys, uint8_t val){
  polyline_t* it = polys;
  float ymin = im->h;
  float ymax = 0;
  while(it){
    point_t* jt = it->head;
    while(jt){
      jt->x += 0.1*(((float)rand()/(float)RAND_MAX)-0.5);
      jt->y += 0.1*(((float)rand()/(float)RAND_MAX)-0.5);
      ymin = fmin(ymin,jt->y);
      ymax = fmax(ymax,jt->y);
      jt = jt->next;
    }
    it->head->x = it->tail->x;
    it->head->y = it->tail->y;
    it = it->next;
  }
  for (int i = ymin; i <= ymax; i++){
    tmp_data_clear();

    polyline_t* it = polys;
    while(it){
      if (fabs(poly_area(it))<4){
        it = it->next;
        continue;
      }
      point_t* jt = it->head;
      while(jt->next){
        point_t p0 = *jt;
        point_t p1;
        // point_t tmp;
        // if (jt->next){
          p1 = *(jt->next);
        // }else{
        //   p1 = *(it->head);
        // }
        point_t pt;
        if (seg_isect_y_line(i,p0,p1,&pt)){
          *(float*)(tmp_data_add(sizeof(float))) = pt.x;
        }
        jt = jt->next;
      }
      it = it->next;
    }
    int npt = tmp_data_len/(sizeof(float));
    qsort(tmp_data, npt, sizeof(float), cmp_ptx);

    // if (i==273){
    //   printf("%d\n",npt);
    //   for (int j = 0; j < npt; j++){
    //     float x =  ( ((float*)tmp_data) [j]);
    //     printf("%f\n",x);
    //   }
    // }

    int j = 0;
    int on = 1;
    while (j < npt-1){
      float x0 =  ( ((float*)tmp_data) [j]);
      float x1 =  ( ((float*)tmp_data) [j+1]);
      if (x0 == x1 && npt % 2){
        if (!on){
          j+=2;
          on = !on;
          continue;
        }else{
          j++;
          continue;
        }
      }
      if (on){
        int x0i = x0;
        int x1i = x1;
        memset(im->data + ((i * im->w )+ x0i), val, x1i-x0i);
      }
      on = !on;
      j ++;
    }
  }
}


float g_rate = 5;
int g_smooth = 5;
int g_resamp = 3;
float g_minarea = 100;
int g_zerosum = 0;

void accentuate_outline(polyline_t* poly, rast_t* rate_mask, int iter){
   
  for (int i = 0; i < iter; i++){   
    if (fabs(poly_area(poly)) < g_minarea){
      break;
    } 
    resample(poly,g_resamp);

    if (poly->size <= 3) break;

    compute_normal(poly);

    compute_curvature(poly);
    smoothen_curvature(poly,g_smooth);

    point_t* it = poly->head;
    if (rate_mask){
      while (it){
        it->curv *= (float)rast_get(rate_mask,it->x,it->y,0)/255.0;
        it = it->next;
      }
    }

    float net = 0;
    it = poly->head;
    while (it){
      net += it->curv;
      it = it->next;
    }
    float err = net/(float)poly->size;

    it = poly->head;
    while (it){
      float c = it->curv;
      if (g_zerosum){
        c -= err;
      }
      it->x += it->nx * c * g_rate;
      it->y += it->ny * c * g_rate;
      it = it->next;
    }
    poly->head->x = poly->tail->x;
    poly->head->y = poly->tail->y;
  }

  // approx_poly_dp(poly, 2);
}

void accentuate_voxels(rast_t* vox, rast_t* rate_mask, int iters, const char* dump_pth){
  static char fname[1024];

  rast_t tmp0 = new_rast(vox->w,vox->h,1);

  for (int z = 0; z < vox->d; z++){

    if (z%50==0) printf("[accentuate] processing layer %d / %d ...\n",z,vox->d);
    rast_t im;
    im.w = vox->w;
    im.h = vox->h;
    im.d = 1;
    im.data = vox->data + ((uint64_t)z*(uint64_t)vox->w*(uint64_t)vox->h);

    rast_t rm;
    if (rate_mask){
      rm.w = rate_mask->w;
      rm.h = rate_mask->h;
      rm.d = 1;
      rm.data = rate_mask->data + (z*rate_mask->w*rate_mask->h);
    }

    for (int i = 0; i < tmp0.w*tmp0.h; i++){
      tmp0.data[i] = !!vox->data[z*vox->w*vox->h+i];
      im.data[i] = 0;
    }
    if (dump_pth){
      sprintf(fname,"%s/%3d_ctr0.gif",dump_pth,z);
      write_gif(fname,&tmp0);
    }
    
    polyline_t* contours = find_contours((int8_t*)tmp0.data, tmp0.w, tmp0.h);
    if (dump_pth){
      sprintf(fname,"%s/%3d_ctr0.svg",dump_pth,z);
      polylines_to_svg(fname,0,0,im.w,im.h, contours);
    }

    polyline_t* it = contours;


    while (it){
      polyline_t* nxt = it->next;
      
      accentuate_outline(it,rate_mask?(&rm):NULL,iters);
      it = nxt;
    }

    draw_polygon_fills(&im, contours, 255);

    if (dump_pth){
      sprintf(fname,"%s/%3d_ctra.svg",dump_pth,z);
      polylines_to_svg(fname,0,0,im.w,im.h, contours);
      sprintf(fname,"%s/%3d_ctra.gif",dump_pth,z);
      write_gif(fname,&im);
    }
    // break;
  }
}

void mush_edits(rast_t* vox0, rast_t* vox, int rad, const char* dump_pth){
  static char fname[1024];

  int ksize = rad*2+1;
  float kern[ksize];

  for (int i = 0; i < vox->w*vox->h*vox->d; i++){
    if (vox->data[i]){
      if (vox0->data[i]){
        vox->data[i] = 127;
      }else{
        vox->data[i] = 255;
      }
    }else{
      if (vox0->data[i]){
        vox->data[i] = 0;
      }else{
        vox->data[i] = 127;
      }
    }
  }

  rast_t voy = new_rast(vox->w,vox->h,vox->d);

  float sigma = 0.3*((float)(ksize-1)*0.5 - 1) + 0.8;
  float ss2 = sigma*sigma*2;
  int i; for (i = 0; i < ksize; i++){
    float x = (float)(i-ksize/2);
    float z = exp(-(x*x)/(ss2))/(2.5066282746*sigma);
    kern[i] = z;
  }

  for (int z = 0; z < vox->d; z++){
    if (z % 100 == 0) printf("[blur] z axis, z = %d / %d\n",z,vox->d);
    for (int y = 0; y < vox->h; y++ ){
      for (int x= 0; x < vox->w; x++ ){
        float sum = 0;
        for (int k = 0; k < ksize; k++ ){
          sum += rast_get_cp_border(vox,x,y,z-rad+k) * kern[k];
        }
        rast_set(&voy,x,y,z,sum);
      }
    }
    if (dump_pth){
      rast_t im;
      im.w = voy.w;
      im.h = voy.h;
      im.d = 1;
      im.data = voy.data + ((uint64_t)z*(uint64_t)voy.w*(uint64_t)voy.h);
      sprintf(fname,"%s/%03d_blur.gif",dump_pth,z);
      write_gif(fname,&im);
    }
  }
  memcpy(vox->data,voy.data,vox->w*vox->h*vox->d);
  free(voy.data);

  for (int i = 0; i < vox->w*vox->h*vox->d; i++){
    if (vox->data[i]>192){
      vox->data[i] = 255;
    }else if (vox->data[i]<64){
      vox->data[i] = 0;
    }else{
      vox->data[i] = vox0->data[i];
    }
  }
}



void print_help(){
  printf("accentuate smaller features on perimeters of a voxel model.                 \n");
  printf("USAGE: accentuate_voxels [options] voxel.bin                                \n");
  printf("OPTIONS:                                                                    \n");
  printf("--output, -o   path/out.bin    output path, voxel binary format:            \n");
  printf("                               w (u32), h (u32), d (u32), data (u8,u8,u8,...\n");
  printf("                               voxel(x,y,z) = data[z*w*h + y*w + x]         \n");
  printf("--iters        5               number of iterations                         \n");
  printf("--rate         5               amount of displacement at each iteration     \n");
  printf("--rate_mask    path/vox.bin    optional: a voxel file giving the rate       \n");
  printf("                               multiplier for each voxel 0->0.0, 255->1.0   \n");
  printf("                               useful for excluding regions from distortion \n");
  printf("--smooth       5               smoothing for curvature calculation          \n");
  printf("--resamp       3               contour resample segment length              \n");
  printf("--blur_z       0               radius of blur on z axis (post-processing)   \n");
  printf("--min_area     100             contours smaller than this will be skipped   \n");
  printf("--zero_sum                     flag: balance outline expand and contract    \n");
  printf("--dump         path/folder     dump intermediate visuals (for debugging)    \n");
  printf("--help                         print this message                           \n");
}


int main(int argc, char** argv){
  srand(0x5EED);

  const char* dump_pth = NULL;
  const char* outp_pth = NULL;
  const char* inpt_pth = NULL;
  const char* rate_pth = NULL;
  int iters = 5;
  int blur_z = 0;

  int i = 1;
  while (i < argc){
    if (!strcmp(argv[i],"--iters")){
      iters = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--rate")){
      g_rate = atof(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--smooth")){
      g_smooth = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--resamp")){
      g_resamp = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--min_area")){
      g_minarea = atof(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--zero_sum")){
      g_zerosum = 1;
      i++;
    }else if (!strcmp(argv[i],"--dump")){
      dump_pth = argv[i+1];
      i+=2;
    }else if (!strcmp(argv[i],"--output") || !strcmp(argv[i],"-o")){
      outp_pth = argv[i+1];
      i+=2;
    }else if (!strcmp(argv[i],"--rate_mask")){
      rate_pth = argv[i+1];
      i+=2;
    }else if (!strcmp(argv[i],"--blur_z")){
      blur_z = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")){
      print_help();
      exit(0);
    }else{
      if (inpt_pth){
        printf("[error] cannot parse commandline argument %s.\n",argv[i]);
        exit(1);
      }else{
        inpt_pth = argv[i];
        i++;
      }
    }
  }
  if (!inpt_pth){
    printf("[warn] no input file. (try '--help' for usage.)\n");
    exit(0);
  }
  

  int pad = ceil(g_rate*(float)iters)+2;

  rast_t vox = read_voxels_from_file_with_padding(inpt_pth, pad, pad);

  
  rast_t rmk;
  if (rate_pth){
    rmk = read_voxels_from_file_with_padding(rate_pth, pad, pad);
  }

  rast_t vox0;
  if (blur_z){
    vox0 = new_rast(vox.w,vox.h,vox.d);
    memcpy(vox0.data,vox.data,vox0.w*vox0.h*vox0.d);
  }

  printf("[voxel] read.\n");
  accentuate_voxels(&vox,rate_pth?(&rmk):NULL,iters,dump_pth);

  if (blur_z){
    mush_edits(&vox0,&vox,blur_z,dump_pth);
  }

  write_voxels_to_file(outp_pth,&vox);
}

