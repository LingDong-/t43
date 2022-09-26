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
  char is_air;
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
void rast_repurpose(rast_t* r, int w, int h, int d){
  if (w*h*d > r->cap){
    r->cap = (uint64_t)w*(uint64_t)h*(uint64_t)d;
    r->data = realloc(r->data,r->cap*sizeof(uint8_t));
  }
  r->w = w;
  r->h = h;
  r->d = d;
}

rast_t segmap;
polyline_t* contour_list[255] = {0};
point_t centroid_list[255] = {0};


void contour_list_to_svg(const char* pth, int w, int h){
  FILE *fp;
  fp = fopen(pth, "w+");
  fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"%d %d %d %d\">\n",w,h,0,0,w,h);

  for (int idx = 0; idx < 255; idx++){
    if (!contour_list[idx]){
      continue;
    }
    polyline_t* it = contour_list[idx];
    int sh = (idx % 3) * 8;
    while(it){
      int color = (rand()%0x7F + 0x7E) << sh;

      fprintf(fp, "<path fill=\"none\" stroke-width=\"1\" stroke=\"#%06x\" stroke-linecap=\"round\" stroke-linejoin=\"round\" d=\"M",color);
      point_t* jt = it->head;
      while(jt){
        fprintf(fp,"%f,%f ",jt->x,jt->y);
        jt = jt->next;
      }
      fprintf(fp, "\" />\n");
      it = it->next;
      // break;
    }
  }
  fprintf(fp, "</svg>\n");
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



#define OPENSET_CAP 4096

int openset[OPENSET_CAP];
int* came_from = NULL;
float* gscore = NULL;
float* fscore = NULL;
int astar_old_size = 0;
int neighbors[4];
int n_neighbors = 0;

polyline_t* reconstruct_path(int w, int current){
  polyline_t* path = new_polyline();
  add_point_to_polyline(path, current%w, current/w);
  while ( came_from[current] != -1){
    current = came_from[current];
    add_point_to_polyline(path, current%w, current/w);
  }
  return path;
}

float heur_dist(int w, int a, int b){
  int x0 = a % w;
  int y0 = a / w;
  int x1 = b % w;
  int y1 = b / w;
  float dx = x1-x0;
  float dy = y1-y0;
  return sqrt(dx*dx+dy*dy);
}

void get_neighbors(rast_t* im, int n){
  int c = n % im->w;
  int r = n / im->w;
  n_neighbors = 0;
  if (c > 0) neighbors[n_neighbors++] = n-1;
  if (c < im->w-1) neighbors[n_neighbors++] = n+1;
  if (r > 0) neighbors[n_neighbors++] = n-im->w;
  if (r < im->h-1) neighbors[n_neighbors++] = n+im->w;
}

void get_strict_neighbors(rast_t* im, int n){
  int c = n % im->w;
  int r = n / im->w;
  n_neighbors = 0;
  if (c > 0 && im->data[n-1]>0) neighbors[n_neighbors++] = n-1;
  if (c < im->w-1 && im->data[n+1]>0) neighbors[n_neighbors++] = n+1;
  if (r > 0 && im->data[n-im->w]>0) neighbors[n_neighbors++] = n-im->w;
  if (r < im->h-1 && im->data[n+im->w]>0) neighbors[n_neighbors++] = n+im->w;
}

polyline_t* astar(rast_t* im, int x0, int y0, int x1, int y1, int stricter){
  
  int w = im->w;
  int h = im->h;
  int w_h = w*h;
  int n_openset = 0;  
  
  if (astar_old_size < w_h){
    came_from = realloc(came_from,w_h*sizeof(int));
    gscore = realloc(gscore,w_h*sizeof(float));
    fscore = realloc(fscore,w_h*sizeof(float));
    astar_old_size = w_h;
  }

  int start = y0*w+x0;
  int goal = y1*w+x1;

  came_from[start] = -1;
  openset[n_openset++] = start;

  
  
  for (int i = 0; i < w_h; i++) gscore[i] = INFINITY;
  for (int i = 0; i < w_h; i++) fscore[i] = INFINITY;
  for (int i = 0; i < w_h; i++) came_from[i] = -1;

  gscore[start] = 0;
  fscore[start] = heur_dist(w,start,goal);
  
  while (n_openset){
    // printf("%d\n",n_openset);
    int amin = -1;
    float fmin = INFINITY;
    for (int i = 0; i < n_openset; ++i){
      float f = fscore[openset[i]];
      if ( f < fmin){
        fmin = f;
        amin = i;
      }
    }
    int current = openset[amin];
    
    if (current == goal){
      return reconstruct_path(w,current);
    }
    n_openset--;
    openset[amin] = openset[n_openset];
    
    if (stricter){
      get_strict_neighbors(im,current);
    }else{
      get_neighbors(im,current);
    }

    for (int i = 0; i < n_neighbors; i++){
      int neighbor = neighbors[i];
      int penalty = 0;
      int val = im->data[neighbor];
      if (val == 1){
        penalty = 900;
      }else if (val == 0){
        penalty = 1000;
      }
      float tentative_gscore = gscore[current] + 1 + penalty;
      // printf("%f %d\n",tentative_gscore,(im->data[neighbor]==1?0:1000));
      if (tentative_gscore < gscore[neighbor]){
        came_from[neighbor] = current;
        gscore[neighbor] = tentative_gscore;
        fscore[neighbor] = gscore[neighbor] + heur_dist(w,neighbor,goal);
        openset[n_openset++] = neighbor;
      }
    }
  }
  return NULL;
}

rast_t astar_tmp;
polyline_t* astar_lowres(rast_t* im, int x0, int y0, int x1, int y1, int stricter){
  int n = 4;
  rast_repurpose(&astar_tmp,im->w/n,im->h/n,1);
  for (int i = 0; i < astar_tmp.h; i++){
    for (int j = 0; j < astar_tmp.w; j++){
      astar_tmp.data[i*astar_tmp.w+j] = im->data[(i*n)*im->w+(j*n)];
    }
  }

  polyline_t* p = astar(&astar_tmp,x0/n,y0/n,x1/n,y1/n,stricter);
  if (!p){
    return NULL;
  }
  approx_poly_dp(p,1.5);
  // if (!p){
  //   char fname[256];
  //   sprintf(fname,"dump/a_%d-%d-%d-%d.gif",x0,y0,x1,y1);
  //   write_gif(fname,&astar_tmp);
  // }
  polyline_t* it = p;
  while (it){
    point_t* jt = it->head;
    while (jt){
      jt->x*=n;
      jt->y*=n;
      jt = jt->next;
    }
    it = it->next;
  }
  return p;
}


polyline_t* retrace_polyline(polyline_t* q){
  polyline_t* p = new_polyline();
  point_t* it;
  for (int k = 0; k < 4; k++){
    reverse_polyline(p);
    it = q->head;
    while(it){
      add_point_to_polyline(p, it->x, it->y);
      it = it->next;
    }
  }
  p->is_air = 1;
  return p;
}


point_t floodfill_q(rast_t* im, int x0, int y0, int black, int white){
  polyline_t* Q = new_polyline();
  add_point_to_polyline(Q,x0,y0);
  rast_set(im,x0,y0,0,white);

  double avg_x = 0;
  double avg_y = 0;
  int cnt = 0;

  while (Q->size > 0){

    point_t* n = Q->head;
    Q->head = n->next;
    Q->size--;

    avg_x += n->x;
    avg_y += n->y;
    cnt++;

    if (rast_get(im,n->x-1,n->y,0)==black){
      add_point_to_polyline(Q, n->x-1,n->y);
      rast_set(im,n->x-1,n->y,0,white);
    }
    if (rast_get(im,n->x+1,n->y,0)==black){
      add_point_to_polyline(Q, n->x+1,n->y);
      rast_set(im,n->x+1,n->y,0,white);
    }
    if (rast_get(im,n->x,n->y-1,0)==black){
      add_point_to_polyline(Q, n->x,n->y-1);
      rast_set(im,n->x,n->y-1,0,white);
    }
    if (rast_get(im,n->x,n->y+1,0)==black){
      add_point_to_polyline(Q, n->x,n->y+1);
      rast_set(im,n->x,n->y+1,0,white);
    }
    free(n);
  }
  destroy_polylines(Q);

  avg_x /= (double)cnt;
  avg_y /= (double)cnt;

  point_t p;
  p.x = avg_x;
  p.y = avg_y;
  return p;
}

void segmentation(rast_t* src, rast_t* dst, int start_idx){
  for (int i = 0; i < src->w*src->h; i++){
    dst->data[i] = src->data[i] ? 254 : 255;
  }
  int idx = start_idx;
  for (int i = 0; i < dst->h; i++){
    for (int j = 0; j < dst->w; j++){
      if (rast_get(dst,j,i,0) == 254){
        point_t p = floodfill_q(dst,j,i,254,idx);
        centroid_list[idx] = p;
        idx++;
      }
    }
  }
}



polyline_t* prepend_polyline_pathfind(rast_t* im, polyline_t* q0, polyline_t* q1){
  // q1 -> p -> q0....
  if (!q0){
    return q1;
  }
  // polyline_t* p = astar(im,q1->tail->x,q1->tail->y,q0->head->x,q0->head->y);

  polyline_t* p = astar_lowres(im,q0->head->x,q0->head->y,q1->tail->x,q1->tail->y, 1);
  if (p){
    // approx_poly_dp(p,1);
    p->is_air = 1;
    p = prepend_polyline(q0,p);
    p = prepend_polyline(p,q1);
    return p;
  }
  // if (stricter){
  //   polyline_t* p = new_polyline();
  //   add_point_to_polyline(p,q1->tail->x,q1->tail->y);
  //   add_point_to_polyline(p,rand()%im->w,rand()%im->h);
  //   add_point_to_polyline(p,q0->head->x,q0->head->y);
  //   p->is_air = 1;
  //   p = prepend_polyline(q0,p);
  //   p = prepend_polyline(p,q1);
  //   return p;
  // }
  p = astar_lowres(im,q0->head->x,q0->head->y,q1->tail->x,q1->tail->y, 0);
  if (p){
    // approx_poly_dp(p,1);
    p->is_air = 1;
    p = prepend_polyline(q0,p);
    p = prepend_polyline(p,q1);
    return p;
  }

  // printf("[warn] pathfinding failure!\n");
  // char fname[256];
  // sprintf(fname,"dump/%f-%f-%f-%f.gif",q0->tail->x,q0->tail->y,q1->head->x,q1->head->y);
  // write_gif(fname,im);
  // exit(0);
  return prepend_polyline(q0,q1);
}


// polyline_t* prepend_polyline_pathfind_from_point(rast_t* im, polyline_t* q, int x0, int y0){
//   polyline_t* p = astar_lowres(im,q->head->x,q->head->y,x0,y0, 1);
//   if (p){
//     // approx_poly_dp(p,1);
//     p->is_air = 1;
//     return prepend_polyline(q,p);
//   }
//   return q;
// }


int step = 1;
int get_contour_idx(polyline_t* c){
  return rast_get(&segmap,c->head->x,c->head->y,0)%255;
}

void crosshatch(rast_t* im, int stride, int dirs){
  polyline_t* seg;
  if (dirs & 1){
    for (int i = 0; i < im->h; i+=stride){
      for (int j = 0; j < im->w+1; j++){
        if (!rast_get(im,j-1,i,0) && rast_get(im,j,i,0)){
          seg = new_polyline();
          add_point_to_polyline(seg,j,i);
        }else if (rast_get(im,j-1,i,0) && !rast_get(im,j,i,0)){
          if (j - seg->head->x < step*2){
            destroy_polylines(seg);
          }else{
            add_point_to_polyline(seg,j-1,i);
            int idx = get_contour_idx(seg);
            contour_list[idx] = prepend_polyline(contour_list[idx], seg);
          }
        }
      }
    }
  }
  if (dirs & 2){
    for (int j = 0; j < im->w; j+=stride){
      for (int i = 0; i < im->h+1; i++){
        if (!rast_get(im,j,i-1,0) && rast_get(im,j,i,0)){
          seg = new_polyline();
          add_point_to_polyline(seg,j,i);
        }else if (rast_get(im,j,i-1,0) && !rast_get(im,j,i,0)){
          if (i - seg->head->y < step*2){
            destroy_polylines(seg);
          }else{
            add_point_to_polyline(seg,j,i-1);
            int idx = get_contour_idx(seg);
            contour_list[idx] = prepend_polyline(contour_list[idx], seg);
          }
        }
      }
    }
  }
}


void crosshatch_pathfind(rast_t* im, rast_t* im_dil, int stride, int dirs){
  polyline_t* seg;
  if (dirs & 1){
    for (int i = 0; i < im->h; i+=stride){
      for (int j = 0; j < im->w+1; j++){
        if (!rast_get(im,j-1,i,0) && rast_get(im,j,i,0)){
          seg = new_polyline();
          add_point_to_polyline(seg,j,i);
        }else if (rast_get(im,j-1,i,0) && !rast_get(im,j,i,0)){
          if (j - seg->head->x < step*2){
            destroy_polylines(seg);
          }else{
            add_point_to_polyline(seg,j-1,i);
            if ((i/stride) % 2){
              reverse_polyline(seg);
            }
            int idx = get_contour_idx(seg);
            contour_list[idx] = prepend_polyline_pathfind(im_dil, contour_list[idx], seg);
          }
        }
      }
    }
  }
  if (dirs & 2){
    for (int j = 0; j < im->w; j+=stride){
      for (int i = 0; i < im->h+1; i++){
        if (!rast_get(im,j,i-1,0) && rast_get(im,j,i,0)){
          seg = new_polyline();
          add_point_to_polyline(seg,j,i);
        }else if (rast_get(im,j,i-1,0) && !rast_get(im,j,i,0)){
          if (i - seg->head->y < step*2){
            destroy_polylines(seg);
          }else{
            add_point_to_polyline(seg,j,i-1);
            if ((j/stride) % 2){
              reverse_polyline(seg);
            }
            int idx = get_contour_idx(seg);
            contour_list[idx] = prepend_polyline_pathfind(im_dil, contour_list[idx], seg);
          }
        }
      }
    }
  }
}


void erode(rast_t* r,rast_t* r1, int k){
  for (int i = 0; i < r->h; i++){
    for (int j = 0; j < r->w; j++){
      rast_set(r1,j,i,0,1);
      for (int ki = -k; ki <= k; ki++){
        if (!rast_get(r,j+ki,i,0)){
          rast_set(r1,j,i,0,0);
          goto nextpix0;
        }
      }
      nextpix0:
      continue;
    }
  }
  for (int i = 0; i < r->h; i++){
    for (int j = 0; j < r->w; j++){
      rast_set(r,j,i,0,1);
      for (int ki = -k; ki <= k; ki++){
        if (!rast_get(r1,j,i+ki,0)){
          rast_set(r,j,i,0,0);
          goto nextpix1;
        }
      }
      nextpix1:
      continue;
    }
  }
}



void dilate(rast_t* r,rast_t* r1, int k){
  for (int i = 0; i < r->h; i++){
    for (int j = 0; j < r->w; j++){
      rast_set(r1,j,i,0,0);
      for (int ki = -k; ki <= k; ki++){
        if (rast_get(r,j+ki,i,0)){
          rast_set(r1,j,i,0,1);
          goto nextpix0;
        }
      }
      nextpix0:
      continue;
    }
  }
  for (int i = 0; i < r->h; i++){
    for (int j = 0; j < r->w; j++){
      rast_set(r,j,i,0,0);
      for (int ki = -k; ki <= k; ki++){
        if (rast_get(r1,j,i+ki,0)){
          rast_set(r,j,i,0,1);
          goto nextpix1;
        }
      }
      nextpix1:
      continue;
    }
  }
}

void dilate_segmentation(rast_t* r,rast_t* r1, int k){
  for (int i = 0; i < r->h; i++){
    for (int j = 0; j < r->w; j++){
      rast_set(r1,j,i,0,255);
      for (int ki = -k; ki <= k; ki++){
        int v = rast_get(r,j+ki,i,0);
        if (v!=255){
          rast_set(r1,j,i,0,v);
          goto nextpix0;
        }
      }
      nextpix0:
      continue;
    }
  }
  for (int i = 0; i < r->h; i++){
    for (int j = 0; j < r->w; j++){
      rast_set(r,j,i,0,255);
      for (int ki = -k; ki <= k; ki++){
        int v = rast_get(r1,j,i+ki,0);
        if (v!=255){
          rast_set(r,j,i,0,v);
          goto nextpix1;
        }
      }
      nextpix1:
      continue;
    }
  }
}


#define MAX_ENTRY 1024

typedef struct {
  char* sec;
  char* key;
  char* val;
} entry_t;

entry_t config[MAX_ENTRY];
int n_config = 0;

void parse_ini(const char* pth){
  FILE* fp = fopen(pth,"r");
  fseek(fp, 0L, SEEK_END);
  int sz = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  char* txt = (char*)malloc(sizeof(char)*sz+1);
  
  int i = 0;
  int c;
  while ((c = fgetc(fp)) != EOF){
    txt[i++] = c == '\r' ? '\n' : c;
  }
  txt[i] = 0;

  int state = 0;
  char* sec = NULL;
  char* key = NULL;
  char* val = NULL;
  int quote = 0;
  int esc = 0;

  for (int i = 0; i < sz; i++){
    if (state == 0){
      if (txt[i] == '['){
        state = 1;
        sec = txt + (i + 1);
      }else if (key == NULL && txt[i] != ' ' && txt[i] != '\n'){
        key = txt + i;
        state = 2;
      }else if (txt[i] == '#'){
        state = 4;
      }
    }else if (state == 1){
      if (txt[i] == ']'){
        txt[i] = 0;
        state = 0;
      }
    }else if (state == 2){
      if (txt[i] == ' '){
        txt[i] = 0;
      }else if (txt[i] == '='){
        int j = i+1;
        do{
          val = txt + j;
          j++;
        }while (val[0] == ' ');
        quote = 0;
        if (val[0] == '"'){
          quote = 1;
          val++;
        }
        state = 3;
      }
    }else if (state == 3){
      if (txt[i] == '\\'){
        esc = 1;
      }else if (!esc){
        if ( ((txt[i] == '\n' || txt[i] == '#') && !quote)  || (txt[i] == '"' && quote && (txt+i+1 != val))){
          if (txt[i] == '#'){
            state = 4;
          }else{
            state = 0;
          }
          txt[i] = 0;
          int j = i;
          while (txt[j] == ' '){
            txt[j--] = 0;
          }
          config[n_config].key = key;
          config[n_config].val = val;
          config[n_config].sec = sec;
          n_config = (n_config+1)%MAX_ENTRY;
          key = NULL;
          val = NULL;

        }
      }else if (esc){
        esc = 0;
      }
    }else if (state == 4){
      if (txt[i] == '\n'){
        state = 0;
      }
    }
  }
}

char* get_config(const char* sec, const char* key){
  for (int i = 0; i < n_config; i++){
    if (sec != config[i].sec && strcmp(sec,config[i].sec)){
      continue;
    }
    if (strcmp(key,config[i].key)){
      continue;
    }
    return config[i].val;
  }
  return NULL;
}

void print_config(){
  for (int i = 0; i < n_config; i++){
    printf("%s.%s=\"%s\"\n",config[i].sec,config[i].key,config[i].val);
  }
}

char* cfg_tmp = NULL;
#define GET_CONFIG_INT(key,defau) int   key = (cfg_tmp = get_config(NULL,#key), cfg_tmp ? atoi(cfg_tmp) : (defau));
#define GET_CONFIG_FLT(key,defau) float key = (cfg_tmp = get_config(NULL,#key), cfg_tmp ? atof(cfg_tmp) : (defau));
#define GET_CONFIG_STR(key,defau) char* key = (cfg_tmp = get_config(NULL,#key), cfg_tmp ?      cfg_tmp  : (defau));
#define GET_CONFIG_CHR(key,defau) char  key = (cfg_tmp = get_config(NULL,#key), cfg_tmp ?    cfg_tmp[0] : (defau));

void rand_rot(polyline_t* poly){
  if (poly->size < 4){
    return;
  }
  point_t* it = poly->head->next;
  point_t* prev = poly->head;
  int n = rand()%(poly->size-2);
  for (int i = 0; i < n; i++){
    prev = prev->next;
    it = it->next;
  }
  poly->tail->next = poly->head;
  poly->head = it;
  poly->tail = prev;
  prev->next = NULL;
}

void reorient_seam(polyline_t* poly, char so){
  if (so == '0'){
    return;
  }
  if (so == 'R'){
    return rand_rot(poly);
  }
  float opt = (so == 'W' || so == 'N') ? INFINITY : -INFINITY;
  point_t* bit = poly->head;
  point_t* bpv = NULL;
  point_t* it = poly->head;
  point_t* pv = NULL;
  while (it){
    float no;
    if (
      (so == 'W' && (no = it->x) < opt) ||
      (so == 'E' && (no = it->x) > opt) ||
      (so == 'N' && (no = it->y) < opt) ||
      (so == 'S' && (no = it->y) > opt)
    ){
      opt = no;
      bit = it;
      bpv = pv;
    }
    pv = it;
    it = it->next;
  }
  if (bpv){
    poly->tail->next = poly->head;
    poly->head = bit;
    poly->tail = bpv;
    bpv->next = NULL;
  }
}




void add_simp_contours(rast_t* r, int dx, int dy, char seam_orientation){
  polyline_t* c = find_contours((int8_t*)r->data, r->w, r->h);
  
  while (c){
    
    polyline_t* c_nxt = c->next;
    c->next = NULL;
    int idx = get_contour_idx(c);

    reorient_seam(c,seam_orientation);
    
    approx_poly_dp(c,0.8);
    if (dx || dy){
      polyline_t* it = c;
      while (it){
        point_t* jt = it->head;
        while (jt){
          jt->x += dx;
          jt->y += dy;
          jt = jt->next;
        }
        it = it->next;
      }
    }
    contour_list[idx] = prepend_polyline(contour_list[idx],c);
    c = c_nxt;
  }
}


void add_simp_contours_pathfind(rast_t* r, rast_t* im_dil, int dx, int dy, char seam_orientation){
  polyline_t* c = find_contours((int8_t*)r->data, r->w, r->h);
  while (c){
    polyline_t* c_nxt = c->next;
    c->next = NULL;

    int idx = get_contour_idx(c);

    reorient_seam(c,seam_orientation);

    approx_poly_dp(c,0.8);
    if (dx || dy){
      polyline_t* it = c;
      while (it){
        point_t* jt = it->head;
        while (jt){
          jt->x += dx;
          jt->y += dy;
          jt = jt->next;
        }
        it = it->next;
      }
    }
    contour_list[idx] = prepend_polyline_pathfind(im_dil,contour_list[idx],c);
    c = c_nxt;
  }
}


void print_header(FILE* fp){
  GET_CONFIG_STR(header,"G90 ; use absolute coordinates\n"
  "M83 ; extruder relative mode\n"
  "M140 S60 ; set final bed temp\n"
  "M104 S150 ; set temporary nozzle temp to prevent oozing during homing\n"
  "G4 S10 ; allow partial nozzle warmup\n"
  "G28 ; home all axis\n"
  "G1 Z50 F240\n"
  "G1 X2 Y10 F3000\n"
  "M104 S220 ; set final nozzle temp\n"
  "M190 S60 ; wait for bed temp to stabilize\n"
  "M109 S220 ; wait for nozzle temp to stabilize\n"
  "G1 Z0.28 F240\n"
  "G92 E0\n"
  "G1 Y140 E10 F1500 ; prime the nozzle\n"
  "G1 X2.3 F5000\n"
  "G92 E0\n"
  "G1 Y10 E10 F1200 ; prime the nozzle\n"
  "G92 E0\n"
  "G21 ; set units to millimeters\n"
  "G90 ; use absolute coordinates\n"
  "M83 ; use relative distances for extrusion\n"
  "; Filament gcode\n"
  "M107\n"
  "G92 E0\n");
  fputs(header,fp);
}

void print_footer(FILE* fp){
  GET_CONFIG_STR(footer,"M107\n"
  "G1 Z200 F600 ; Move print head up\n"
  "G1 X5 Y182.4 F9000 ; present print\n"
  "G1 Z220 F600 ; Move print head further up\n"
  "G1 Z200 F600 ; Move print head further up\n"
  "M140 S0 ; turn off heatbed\n"
  "M104 S0 ; turn off temperature\n"
  "M107 ; turn off fan\n"
  "M84 X Y E ; disable motors\n");
  fputs(footer,fp);
}

void print_layer_start(FILE* fp, float z, int feedrate){
  fprintf(fp, "G1 Z%f F%d\n",z,feedrate);
}

void print_layer_end(FILE* fp){
  fprintf(fp, "G92 E0\n");
}

inline float dist2d(float x0, float y0, float x1, float y1){
  float dx = x1-x0;
  float dy = y1-y0;
  return sqrt(dx*dx+dy*dy);
}


void voxels_to_gcode(FILE* fp, rast_t* r, const char* dump_pth){

  static char fname[1024];

  GET_CONFIG_INT(num_shells,3);
  GET_CONFIG_INT(num_brims,12);
  GET_CONFIG_INT(num_ceilings,3);
  GET_CONFIG_INT(support_separation,7);
  GET_CONFIG_INT(infill,16);
  GET_CONFIG_INT(support_fill,8);
  GET_CONFIG_FLT(layer_height,0.2);
  GET_CONFIG_FLT(extrusion_factor,0.032);
  GET_CONFIG_FLT(extrusion_factor_lay0,0.028);
  GET_CONFIG_FLT(origin_x,50);
  GET_CONFIG_FLT(origin_y,50);

  GET_CONFIG_INT(feedrate,900);
  GET_CONFIG_INT(feedrate_lay0,1600);

  GET_CONFIG_CHR(seam_orientation,'0');

  GET_CONFIG_FLT(retract,-3);
  GET_CONFIG_FLT(retract_more,-1);
  GET_CONFIG_FLT(unretract,5);
  GET_CONFIG_FLT(retract_speed,2);
  GET_CONFIG_INT(min_jump,3);
  GET_CONFIG_FLT(min_wipe,2);


  GET_CONFIG_FLT(shell_offset_factor,1.5);
  GET_CONFIG_FLT(shell0_offset_factor,0.5);
  

  support_separation *= step;
  infill *= step;
  support_fill *= step;

  int shell_ofs = (int)ceil((float)step*shell_offset_factor);

  rast_t tmp0 = new_rast(r->w,r->h,1);
  rast_t tmp1 = new_rast(r->w,r->h,1);
  rast_t tmp2 = new_rast(r->w,r->h,1);
  int brimpad = num_brims*shell_ofs+1;
  rast_t tmp3 = new_rast(r->w+brimpad*2,r->h+brimpad*2,1);
  rast_t tmp4 = new_rast(r->w+brimpad*2,r->h+brimpad*2,1);
  rast_t tmp5 = new_rast(r->w,r->h,1);
  rast_t tmp6 = new_rast(r->w,r->h,1);
  rast_t tmp7 = new_rast(r->w,r->h,1);
  rast_t tmp8 = new_rast(r->w,r->h,1);
  
  rast_repurpose(&segmap,r->w,r->h,1);

  print_header(fp);

  float prvlay_x;
  float prvlay_y;

  for (int z = 0; z < r->d; z+=step){
    if (z % (step*10) == 0){
      printf("[slicing] %d / %d\n",z,r->d);
    }
    rast_t im;
    im.w = r->w;
    im.h = r->h;
    im.d = 1;
    im.data = r->data + ((uint64_t)z*(uint64_t)r->w*(uint64_t)r->h);

    for (int i = 0; i < r->w*r->h; i++){
      tmp0.data[i] = im.data[i]==255;
      tmp5.data[i] = im.data[i]==255;
      tmp6.data[i] = im.data[i]==255;
      tmp7.data[i] = im.data[i]==255;
    }

    // segmentation(&im,&segmap);
    // dilate_segmentation(&segmap,&tmp1,4);

    for (int i = 0; i < min_jump; i++){
      dilate(&tmp7,&tmp1,step);
    }

    segmentation(&tmp7,&segmap,0);
    dilate_segmentation(&segmap,&tmp1,1);

    // PREP SUPPORT -------------------
    if (support_fill){
      for (int i = 0; i < r->w*r->h; i++){
        tmp2.data[i] = im.data[i] == 255;
      }
      dilate(&tmp2,&tmp1,support_separation);
      for (int i = 0; i < r->w*r->h; i++){
        tmp1.data[i] = im.data[i]==2 && !tmp2.data[i];
      }
      segmentation(&tmp1,&tmp2,128);
      for (int i = 0; i < r->w*r->h; i++){
        if (tmp2.data[i] != 255 && segmap.data[i] == 255){
          segmap.data[i] = tmp2.data[i];
        }
      }
    }

    if (dump_pth){
      sprintf(fname,"%s/%3d_semt.gif",dump_pth,z);
      write_gif(fname,&segmap);
    }
  
    erode(&tmp5,&tmp1,shell_ofs);
    dilate(&tmp6,&tmp1,1);
    for (int i = 0; i < r->w*r->h; i++){
      tmp5.data[i] += tmp6.data[i];
    }

    // SHELLS -------------------------
    if (z > 0 || num_brims == 0){
      for (int k = 0; k < num_shells; k++){
        if (k){
          erode(&tmp0,&tmp1,shell_ofs);
        }else{
          erode(&tmp0,&tmp1,(int)((float)step*shell0_offset_factor));
        }
        add_simp_contours_pathfind(&tmp0,&tmp5,0,0,seam_orientation);
      }
    }

    // CEILINGS -------------------------

    rast_clear(&tmp2);
    for (int i = 0; i < r->h; i++){
      for (int j = 0; j < r->w; j++){
        for (int k = 1; k < num_ceilings+1; k++){
          if (rast_get(&tmp0,j,i,0)==1 && rast_get(r,j,i,z+step*k)!=255){
            rast_set(&tmp2,j,i,0,1);
          }
          if (rast_get(&tmp0,j,i,0)==1 && rast_get(r,j,i,z-step*k)!=255){
            rast_set(&tmp2,j,i,0,1);
          }
        }
      }
    }
    // erode(&tmp2,&tmp1,2);

    // contours = crosshatch(&tmp2,contours,shell_ofs*(z?1:2),1+(z/step)%2);
    if (z){
      crosshatch_pathfind(&tmp2,&tmp5,shell_ofs,1+(z/step)%2);
    }else{
      crosshatch(&tmp2,shell_ofs*(z?1:2),1+(z/step)%2);
    }
    // INFILL -------------------------

    for (int i = 0; i < r->w*r->h; i++){
      if (tmp0.data[i] && tmp2.data[i]){
        tmp0.data[i] = 0;
      }
    }
    // if (z == 0){
      dilate(&tmp0,&tmp1,1);
    // }

    
    crosshatch_pathfind(&tmp0,&tmp5,infill,3);

    // SUPPORT -------------------------

    // erode(&tmp0,&tmp1,support_separation);
    // contours = add_simp_contours(&tmp0,contours,0,0,seam_orientation);
    // contours = crosshatch(&tmp0,contours,support_fill,2-(z/step)%2);

    if (support_fill){
      for (int i = 0; i < r->w*r->h; i++){
        if (segmap.data[i]>=128 && segmap.data[i] != 255 && r->data[MIN(z+step,r->d-1)*r->w*r->h+i] != 255){
          tmp0.data[i] = 1;
          tmp2.data[i] = 1;
        }else{
          tmp0.data[i] = 0;
          tmp2.data[i] = 0;
        }
      }
      dilate(&tmp2,&tmp1,1);
      crosshatch_pathfind(&tmp0,&tmp2,support_fill,(z/step)%3);
    }

    // BRIM -------------------------

    if (z == 0){
      for (int i = 0; i < r->h; i++){
        for (int j = 0; j < r->w; j++){
          rast_set(&tmp3,j+brimpad,i+brimpad,0,!!rast_get(&im,j,i,0));
        }
      }
      for (int k = 0; k < num_brims+1; k++){

        if (k){
          dilate(&tmp3,&tmp4,shell_ofs);
        }else{
          dilate(&tmp3,&tmp4,1);
          erode(&tmp3,&tmp4,1);
        }
        // write_gif("brim.gif",&tmp3);
        add_simp_contours(&tmp3,-brimpad,-brimpad,seam_orientation);
      }
    
    }

    for (int i = 0; i < r->w*r->h; i++){
      int v = MAX(tmp7.data[i],tmp5.data[i]);
      tmp7.data[i] = MAX(v, tmp8.data[i]);
      tmp8.data[i] = v;
    }


    print_layer_start(fp, (z/step)*layer_height+layer_height, ((z/step) < num_ceilings) ? feedrate_lay0 : feedrate );


    float x0;
    float y0;

    float lcx;
    float lcy;

    for (int idx = 0; idx < 255; idx++){
      if (!contour_list[idx]){
        continue;
      }
      
      int tryfind = !idx;
      if (tryfind && z){
        polyline_t* p = astar_lowres(&tmp7,contour_list[idx]->head->x,contour_list[idx]->head->y,prvlay_x,prvlay_y, 1);
        if (p){
          p->is_air = 1;
          contour_list[idx] = prepend_polyline(contour_list[idx],p);
        }else{
          tryfind = 0;
        }
      }

      polyline_t* it = contour_list[idx];
  
      while (it){

        point_t* jt = it->head;
        while (jt){
          float x = origin_x+(jt->x)/(float)step*layer_height;
          float y = origin_y+(jt->y)/(float)step*layer_height;

          if (jt == it->head && it == contour_list[idx]){
            if (!tryfind){
              // printf("forced jump at %d %f\n",z,(z/step)*layer_height+layer_height);
              if (dist2d(lcx,lcy,x0,y0) > min_wipe){
                fprintf(fp, "G1 X%f Y%f E%f F%d\n",lcx,lcy,retract,(int)(feedrate*retract_speed));
              }else{
                
                float rx = origin_x+(rand()%im.w)/(float)step*layer_height;
                float ry = origin_y+(rand()%im.h)/(float)step*layer_height;
                
                fprintf(fp, "G1 X%f Y%f E%f F%d\n",rx,ry,retract,(int)(feedrate*retract_speed));
              }
              
              fprintf(fp, "G1 X%f Y%f E%f F%d\n",x,y,retract_more,feedrate);
              fprintf(fp, "G1 E%f F%d\n",unretract,(int)(feedrate*retract_speed));
              fprintf(fp, "G1 F%d\n",feedrate);
            }
          }else{
            float d = dist2d(x,y,x0,y0);
            fprintf(fp, "G1 X%f Y%f E%f\n",x,y,d*(z?extrusion_factor:extrusion_factor_lay0));
          }

          prvlay_x = jt->x;
          prvlay_y = jt->y;
          x0 = x;
          y0 = y;
          jt = jt->next;
        }
        
        it = it->next;
      }

      lcx = origin_x+(centroid_list[idx].x)/(float)step*layer_height;
      lcy = origin_y+(centroid_list[idx].y)/(float)step*layer_height;

    }

    if (!z){
      fprintf(fp,"M106 S255\n");
    }
    print_layer_end(fp);

    if (dump_pth != NULL){
      sprintf(fname,"%s/%3d_cntr.svg",dump_pth,z);
      contour_list_to_svg(fname,r->w,r->h);
    }

    for (int idx = 0; idx < 255; idx++){
      if (contour_list[idx]){
        destroy_polylines(contour_list[idx]);
        contour_list[idx] = NULL;
      }
    }
  }
  print_footer(fp);

  rast_destroy(&tmp0);
  rast_destroy(&tmp1);
  rast_destroy(&tmp2);
  rast_destroy(&tmp3);
  rast_destroy(&tmp4);
  rast_destroy(&tmp5);
  rast_destroy(&tmp6);
  rast_destroy(&tmp7);
  rast_destroy(&tmp8);
}


void print_help(){
  printf("voxel to gcode, for 3d printing. 255=solid, 2=support                       \n");
  printf("USAGE: voxels_to_gcode [options] input.bin                                  \n");
  printf("voxel binary format: w (u32 LE), h (u32), d (u32), data (u8,u8,u8,...       \n");
  printf("voxel(x,y,z) = data[z*w*h + y*w + x]    0=void, 255=solid, 2=support        \n");
  printf("OPTIONS:                                                                    \n");
  printf("--output, -o   path/out.gcode  output path                                  \n");
  printf("--config       path/cfg.ini    configuration file for all parameters        \n");
  printf("                               (see profiles/ender3.ini for example)        \n");
  printf("--step         1               step = n  ->  layer[i] = voxel(:,:,n*i)      \n");
  printf("                               downscale voxel, but give more XY resolution \n");
  printf("                               recommended: 4                               \n");
  printf("--sink         0               number of layers to sink into build plane (z)\n");
  printf("--dump         path/folder     dump intermediate visuals (for debugging)    \n");
  printf("--help                         print this message                           \n");
}

int main(int argc, char** argv){
  srand(0x5EED);

  const char* dump_pth = NULL;
  const char* outp_pth = NULL;
  const char* inpt_pth = NULL;
  const char* cnfg_pth = NULL;
  int sink = 0;

  int i = 1;
  while (i < argc){
    if (!strcmp(argv[i],"--config")){
      cnfg_pth = argv[i+1];
      i+=2;
    }else if (!strcmp(argv[i],"--output") || !strcmp(argv[i],"-o")){
      outp_pth = argv[i+1];
      i+=2;
    }else if (!strcmp(argv[i],"--dump") ){
      dump_pth = argv[i+1];
      i+=2;
    }else if (!strcmp(argv[i],"--step")){
      step = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--sink")){
      sink = atoi(argv[i+1]);
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
  if (cnfg_pth == NULL){
    cnfg_pth = "profiles/ender3.ini";
  }

  parse_ini(cnfg_pth);
  // print_config();
  rast_t vox = read_voxels_from_file(inpt_pth);
  printf("[voxel] loaded.\n");

  rast_t vox2;
  vox2.w = vox.w;
  vox2.h = vox.h;
  vox2.d = vox.d-sink;
  vox2.data = vox.data + ((uint64_t)sink*(uint64_t)vox.w*(uint64_t)vox.h);

  FILE* fp = fopen(outp_pth,"w+");
  voxels_to_gcode(fp,&vox2,dump_pth);
  fclose(fp);
  return 0;
}
