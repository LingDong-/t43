#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct {
  float x;
  float y;
  float z;
}vec3_t;

typedef struct {
  vec3_t* vtxs;
  int* tris;
  int n_vtx;
  int n_tri;
}mesh_t;

typedef struct {
  vec3_t min;
  vec3_t max;
}bbox_t;

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
  int lvl;
} polyline_t;

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
uint8_t rast_get(rast_t* r, int x, int y, int z){
  if (x < 0 || x >= r->w || y < 0 || y >= r->h || z < 0 || z >= r->d){
    return 0;
  }
  return r->data[z*(r->w*r->h)+y*(r->w)+x];
}
uint8_t rast_get_white_border(rast_t* r, int x, int y, int z){
  if (x < 0 || x >= r->w || y < 0 || y >= r->h || z < 0 || z >= r->d){
    return 255;
  }
  return r->data[z*(r->w*r->h)+y*(r->w)+x];
}
uint8_t rast_set(rast_t* r, int x, int y, int z, uint8_t v){
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
void rast_repurpose(rast_t* r, int w, int h, int d){
  if (w*h*d > r->cap){
    r->cap = (uint64_t)w*(uint64_t)h*(uint64_t)d;
    r->data = realloc(r->data,r->cap*sizeof(uint8_t));
  }
  r->w = w;
  r->h = h;
  r->d = d;
}

#define BLEND_MIN 1
#define BLEND_MAX 2

void rast_blend(rast_t* r0, rast_t* r1, int mode){
  for (int i = 0; i < r0->w*r0->h*r0->d; i++){
    if (mode == BLEND_MAX){
      r0->data[i] = MAX(r0->data[i],r1->data[i]);
    }else if (mode == BLEND_MIN){
      r0->data[i] = MIN(r0->data[i],r1->data[i]);
    }
  }
}


polyline_t* new_polyline(){
  polyline_t* q0 = (polyline_t*)malloc(sizeof(polyline_t));
  q0->head = NULL;
  q0->tail = NULL;
  q0->prev = NULL;
  q0->next = NULL;
  q0->size = 0;
  q0->lvl = -1;
  return q0;
}

void print_polyline(polyline_t* q){
  if (!q){
    return;
  }
  point_t* jt = q->head;
  while(jt){
    printf("%f,%f ",jt->x,jt->y);
    jt = jt->next;
  }
  printf("\n");
}

void polylines_to_svg(const char* pth, int xmin, int ymin, int xmax, int ymax, polyline_t* q){
  FILE *fp;
  fp = fopen(pth, "w+");
  int w = xmax-xmin;
  int h = ymax-ymin;
  fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"%d %d %d %d\">\n",w*10,h*10,xmin,ymin,w,h);

  polyline_t* it = q;
  while(it){
    fprintf(fp, "<path fill=\"none\" stroke-width=\"0.1\" stroke=\"#%06x\" stroke-linecap=\"round\" stroke-linejoin=\"round\" d=\"M",rand()%0xFFFFFF);
    point_t* jt = it->head;
    while(jt){
      fprintf(fp,"%f,%f ",jt->x,jt->y);
      jt = jt->next;
    }
    fprintf(fp, "z\" />\n");
    it = it->next;
    // break;
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
    fputc(i*2,fp);fputc(i*2,fp);fputc(i*2,fp);
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
      fputc(im->data[i*126+j]/2,fp);
    }
  }
  if (exc){
    fputc(exc+1,fp);
    fputc(0x80,fp);
    for (int i = 0; i < exc; i++){
      fputc(im->data[inc+i]/2,fp);
    }
  }
  fputc(0x01,fp);fputc(0x81,fp);fputc(0x00,fp);fputc(0x3B,fp);
  fclose(fp);
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

void cat_head_polyline(polyline_t* q0, polyline_t* q1){
  if (!q1){
    return;
  }
  if (!q0){
    *q0 = *new_polyline();
  }
  if (!q1->head){
    return;
  }
  if (!q0->head){
    q0->head = q1->head;
    q0->tail = q1->tail;
    return;
  }
  q1->tail->next = q0->head;
  q0->head = q1->head;
  q0->size += q1->size;
  q0->tail->next = NULL;
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

mesh_t load_obj(const char* path){
  mesh_t obj;

  FILE * fp;
  char * line = NULL;
  size_t len;
  size_t linelen = 0;
  size_t read;
  fp = fopen(path,"r");

  int nv = 0;
  int nf = 0;
  while ((read = getline(&line, &len, fp)) != -1) {
    linelen = strlen(line)-1;
    if (linelen < 2){
      continue;
    }
    if (line[0] == 'v' && line[1] == ' '){
      nv ++;
    }else if (line[0] == 'f' && line[1] == ' '){
      nf ++;
    }
  }

  obj.n_tri = nf;
  obj.n_vtx = nv;

  obj.vtxs = (vec3_t*)malloc(sizeof(vec3_t)*nv);

  obj.tris = (int*)malloc(sizeof(int)*nf*3);

  rewind(fp);

  linelen = 0;

  int vi = 0;
  int fi = 0;

  while ((read = getline(&line, &len, fp)) != -1) {
    linelen = strlen(line)-1;
    if (linelen < 2){
      continue;
    }
    if (line[0] == 'v' && line[1] == ' '){
      float x, y, z /*-Wall*/ = 0.0 /**/;
      int xyz = 0;
      int xi = 2; int yi = 2; int zi = 2;
      int i; for (i= 2; i < linelen+1; i++){
        if (line[i] == ' ' || i == linelen){
          line[i] = 0;
          if (xyz == 0){
            x = atof(&line[xi]);
            yi = i+1;
          }else if (xyz == 1){
            y = atof(&line[yi]);
            zi = i+1;
          }else if (xyz == 2){
            z = atof(&line[zi]);
            break;
          }
          xyz ++;
        }
      }
      obj.vtxs[vi].x = x;
      obj.vtxs[vi].y = y;
      obj.vtxs[vi].z = z;
      // printf("%f %f %f %d %d %d\n",x,y,z,xi,yi,zi);
      vi ++;
    }else if (line[0] == 'f' && line[1] == ' '){
      int a,b,c /*-Wall*/ = 0.0 /**/;
      int abc = 0;
      int ai = 2; int bi = 2; int ci = 2;
      int i; for (i= 2; i < linelen+1; i++){
        if (line[i] == '/'){
          line[i] = 0;
        }if (line[i] == ' ' || i == linelen){
          line[i] = 0;
          if (abc == 0){
            a = atoi(&line[ai]);
            bi = i+1;
          }else if (abc == 1){
            b = atoi(&line[bi]);
            ci = i+1;
          }else if (abc == 2){
            c = atoi(&line[ci]);
            break;
          }
          abc ++;
        }
      }
      a += (a < 0) ? vi : (-1);
      b += (b < 0) ? vi : (-1);
      c += (c < 0) ? vi : (-1);
      obj.tris[fi*3  ] = a;
      obj.tris[fi*3+1] = b;
      obj.tris[fi*3+2] = c;
      // printf("%d %d %d\n",a,b,c);
      fi++;
    }
  }
  return obj;
}


bbox_t bounding_box(mesh_t* mesh){
  bbox_t bb;
  bb.min.x = INFINITY;
  bb.min.y = INFINITY;
  bb.min.z = INFINITY;
  bb.max.x = -INFINITY;
  bb.max.y = -INFINITY;
  bb.max.z = -INFINITY;
  for (int i = 0; i < mesh->n_vtx; i++){
    bb.min.x = fmin(bb.min.x, mesh->vtxs[i].x);
    bb.min.y = fmin(bb.min.y, mesh->vtxs[i].y);
    bb.min.z = fmin(bb.min.z, mesh->vtxs[i].z);
    bb.max.x = fmax(bb.max.x, mesh->vtxs[i].x);
    bb.max.y = fmax(bb.max.y, mesh->vtxs[i].y);
    bb.max.z = fmax(bb.max.z, mesh->vtxs[i].z);
  }
  return bb;
}

float dist2d(point_t* a, point_t *b){
  float dx = a->x - b->x;
  float dy = a->y - b->y;
  return sqrt(dx*dx+dy*dy);
}

void pop_head_polyline(polyline_t* q){
  point_t* p = q->head;
  q->head = q->head->next;
  free(p);
  q->size --;
}

polyline_t* connect_segs(polyline_t* segs,float eps){
  int flag = 0;
  polyline_t* it = segs;

  while(it){
    if (it->size == 0){
      it = it->next;
      continue;
    }
    point_t* p0 = it->head;
    point_t* p1 = it->tail;

    polyline_t* jt = it->next;
    int cnt[4] = {0,0,0,0};
    polyline_t* bjt[4] = {NULL,NULL,NULL,NULL};
    int typ;

    while (jt){
      if (it == jt){
        jt = jt->next;
        continue;
      }
      if (jt->size == 0){
        jt = jt->next;
        continue;
      }
      point_t* q0 = jt->head;
      point_t* q1 = jt->tail;
      if (dist2d(p1,q0)<=eps){
        bjt[0] = jt;
        cnt[0] ++;
      }else if (dist2d(p1,q1)<=eps){
        bjt[1] = jt;
        cnt[1] ++;
      }else if (dist2d(p0,q0)<=eps){
        bjt[2] = jt;
        cnt[2] ++;
      }else if (dist2d(p0,q1)<=eps){
        bjt[3] = jt;
        cnt[3] ++;
      }
      jt = jt->next;
    }
    int ok = 0;
    for (int i = 0; i < 4; i++){
      if (cnt[i] == 1){
        jt = bjt[i];
        typ = i;
        ok = 1;
      }
    }
    
    if (!ok){
      it = it->next;
      continue;
    }
    // printf("%d\n",typ);
    if (typ == 0){
      pop_head_polyline(jt);
      cat_tail_polyline(it, jt);
      jt->size = 0;
      flag = 1;
    }else if (typ == 1){
      reverse_polyline(jt);
      pop_head_polyline(jt);
      cat_tail_polyline(it, jt);

      jt->size = 0;
      flag = 1;
    }else if (typ == 2){
      // reverse_polyline(it);
      // pop_head_polyline(jt);
      // cat_tail_polyline(it,jt);

      reverse_polyline(jt);
      pop_head_polyline(it);
      cat_head_polyline(it,jt);
      jt->size = 0;
      flag = 1;

    }else if (typ == 3){
      pop_head_polyline(it);
      cat_head_polyline(it,jt);
      jt->size = 0;
      flag = 1;
    }
    it = it->next;
  }
  if (!flag){
    return segs;
  }
  it = segs;
  while (it){
    if (it->size == 0){
      if (it->prev){
        it->prev->next = it->next;
      }else{
        segs = it->next;
      }
      if (it->next) it->next->prev = it->prev;
      free(it);
    }
    it = it->next;
  }
  return connect_segs(segs,eps);
}

int seg_isect_z_plane(float z,vec3_t p0,vec3_t p1, vec3_t* out){
  if (p0.z > p1.z){
    vec3_t tmp;
    tmp = p0;
    p0 = p1;
    p1 = tmp;
  }
  if (p1.z == p0.z){
    p1.z+=0.00001;
  }
  float t = (z - p0.z)/(p1.z-p0.z);
  if (t < 0 || t >= 1){
    return 0;
  }
  out->x = p0.x * (1-t) + p1.x * t;
  out->y = p0.y * (1-t) + p1.y * t;
  out->z = z;
  return 1;
}


int line_isect(float p0x, float p0y, float p1x, float p1y, float q0x, float q0y, float q1x, float q1y, float* t, float* s) {
  float d0x = p1x - p0x;
  float d0y = p1y - p0y;
  float d1x = q1x - q0x;
  float d1y = q1y - q0y;
  float vc = d0x * d1y - d0y * d1x;
  if (vc == 0) {
    return 0;
  }
  float vcn = vc * vc;
  float q0x_p0x = q0x - p0x;
  float q0y_p0y = q0y - p0y;
  float vc_vcn = vc / vcn;
  (*t) = (q0x_p0x * d1y - q0y_p0y * d1x) * vc_vcn;
  (*s) = (q0x_p0x * d0y - q0y_p0y * d0x) * vc_vcn;
  return 1;
}

int pt_in_poly(point_t* p,polyline_t* poly){
  int n = 0;
  float qx = p->x + M_PI;
  float qy = p->y + M_E;
  point_t* it = poly->head;
  while (it){
    float x0 = it->x;
    float y0 = it->y;
    float x1, y1;
    if (it->next){
      x1 = it->next->x;
      y1 = it->next->y;
    }else{
      x1 = poly->head->x;
      y1 = poly->head->y;
    }
    float t;
    float s;
    if (line_isect(p->x,p->y,qx,qy,x0,y0,x1,y1,&t,&s)){
      if (t >= 0 && 0 <= s && s <= 1){
        n++;
      }
    }
    it = it->next;
  }
  return n % 2 == 1;
}


int poly_check_internal(polyline_t* polys, polyline_t* p){
  if (p->lvl >= 0){
    return p->lvl;
  }
  polyline_t* it = polys;
  while (it){
    if (it == p){
      it = it->next;
      continue;
    }
    if (pt_in_poly(p->head,it)){
      poly_check_internal(polys,it);
      return (p->lvl = it->lvl+1);
    }
    it = it->next;
  }
  return (p->lvl = 0);
}

int polys_check_internal(polyline_t* polys){
  polyline_t* it = polys;
  int mx = 0;
  while (it){
    mx = MAX(mx,poly_check_internal(polys,it));
    it = it->next;
  }
  return mx;
}



void plot_line(rast_t* im, int x0, int y0, int x1, int y1, int thick, uint8_t val){
  int dx = abs(x1 - x0);
  int sx = x0 < x1 ? 1 : -1;
  int dy = -abs(y1 - y0);
  int sy = y0 < y1 ? 1 : -1;
  int error = dx + dy;
  int ws[thick*2+1];
  for (int i = 0; i < thick*2+1; i++){
    float t = (float)i/(float)(thick*2) * 2 - 1;
    float y = sqrt(1-t*t);
    ws[i] = ceil(y*thick);
  }
  while (1){
    if (!thick){
      rast_set(im, x0, y0, 0, val);
    }else{
      int yl = MAX(y0-thick,0);
      int yr = MIN(y0+thick+1,im->h);
      for (int i = yl; i < yr; i++){
        int ww = ws[i-y0+thick];
        int xl = MAX(x0-ww,0);
        int xr = MIN(x0+1+ww,im->w);
        memset(im->data + (i * im->w) + xl, val, xr-xl);
      }
    }
    if (x0 == x1 && y0 == y1) break;
    int e2 = 2 * error;
    if (e2 >= dy){
      if (x0 == x1) break;
      error += dy;
      x0 += sx;
    }
    if (e2 <= dx){
      if (y0 == y1) break;
      error += dx;
      y0 += sy;
    }
  }
}


#define w2c_x(X) (((X)-bb.min.x)/(bb.max.x-bb.min.x)*(float)im->w)
#define w2c_y(Y) (((Y)-bb.min.y)/(bb.max.y-bb.min.y)*(float)im->h)

void draw_polygon_outlines(rast_t* im, polyline_t* polys, bbox_t bb, int lvl, int thick, uint8_t val){
  polyline_t* it = polys;
  while(it){
    if (it->lvl != lvl && lvl != -1){
      it = it->next;
      continue;
    }
    point_t* jt = it->head;
    while(jt){
      float x0 = (jt->x);
      float y0 = (jt->y);
      float x1,y1;
      if (jt->next){
        x1 = (jt->next->x);
        y1 = (jt->next->y);
      }else{
        x1 = (it->head->x);
        y1 = (it->head->y);
      }
      plot_line(im,round(w2c_x(x0)),round(w2c_y(y0)),round(w2c_x(x1)),round(w2c_y(y1)),thick,val);
      jt = jt->next;
    }
    it = it->next;
  }
}

int seg_isect_y_line(float y, point_t p0, point_t p1, point_t* out){
  if (p0.y > p1.y){
    point_t tmp = p0;
    p0 = p1;
    p1 = tmp;
  }
  if (p0.y == p1.y){
    p1.y += 0.00001;
  }
  float t = (y - p0.y)/(p1.y-p0.y);
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

void draw_polygon_fills(rast_t* im, polyline_t* polys, bbox_t bb, int lvl, uint8_t val){
  polyline_t* it = polys;
  float yminf = INFINITY;
  float ymaxf = -INFINITY;
  while(it){
    if (it->lvl != lvl && lvl != -1){
      it = it->next;
      continue;
    }
    point_t* jt = it->head;
    while(jt){
      yminf = fmin(yminf,jt->y);
      ymaxf = fmax(ymaxf,jt->y);
      jt = jt->next;
    }
    it = it->next;
  }
  int ymin = w2c_y(yminf);
  int ymax = w2c_y(ymaxf);
  // printf("%d %d\n",ymin,ymax);
  // int ymin = 0;
  // int ymax = im->h-1;
  for (int i = ymin; i <= ymax; i++){
    tmp_data_clear();

    polyline_t* it = polys;
    while(it){
      if (it->lvl != lvl && lvl != -1){
        it = it->next;
        continue;
      }
      point_t* jt = it->head;
      while(jt){
        point_t p0 = *jt;
        point_t p1;
        // point_t tmp;
        if (jt->next){
          p1 = *(jt->next);
        }else{
          p1 = *(it->head);
        }
        p0.x = w2c_x(p0.x);
        p0.y = w2c_y(p0.y);
        p1.x = w2c_x(p1.x);
        p1.y = w2c_y(p1.y);
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

    // if (i==201){
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
        // printf("damn!!\n");
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


int thinning_zs_iteration(rast_t* im, int iter) {
  // http://agcggs680.pbworks.com/f/Zhan-Suen_algorithm.pdf
  int W = im->w;
  int H = im->h;
  int diff = 0;
  for (int i = 1; i < H-1; i++){
    for (int j = 1; j < W-1; j++){
      int p2 = im->data[(i-1)*W+j]   & 1;
      int p3 = im->data[(i-1)*W+j+1] & 1;
      int p4 = im->data[(i)*W+j+1]   & 1;
      int p5 = im->data[(i+1)*W+j+1] & 1;
      int p6 = im->data[(i+1)*W+j]   & 1;
      int p7 = im->data[(i+1)*W+j-1] & 1;
      int p8 = im->data[(i)*W+j-1]   & 1;
      int p9 = im->data[(i-1)*W+j-1] & 1;
      
      int A  = (p2 == 0 && p3 == 1) + (p3 == 0 && p4 == 1) +
        (p4 == 0 && p5 == 1) + (p5 == 0 && p6 == 1) +
        (p6 == 0 && p7 == 1) + (p7 == 0 && p8 == 1) +
        (p8 == 0 && p9 == 1) + (p9 == 0 && p2 == 1);
      int B  = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
      int m1 = iter == 0 ? (p2 * p4 * p6) : (p2 * p4 * p8);
      int m2 = iter == 0 ? (p4 * p6 * p8) : (p2 * p6 * p8);
      if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)
        im->data[i*W+j] |= 2;
    }
  }
  for (int i = 0; i < H*W; i++){
    int marker = im->data[i]>>1;
    int old = im->data[i]&1;
    im->data[i] = old & (!marker);
    if ((!diff) && (im->data[i] != old)){
      diff = 1;
    }
  }
  return diff;
};

void thinning_zs(rast_t* im){
  int diff = 1;
  do {
    diff &= thinning_zs_iteration(im,0);
    diff &= thinning_zs_iteration(im,1);
  }while (diff);
}


rast_t tmp_im;
void draw_skeleton(rast_t* im, int thick, uint8_t val){
  int w = im->w/2;
  int h = im->h/2;
  if (!tmp_im.data){
    tmp_im = new_rast(w,h,1);
  }
  rast_repurpose(&tmp_im,w,h,1);
  for (int i = 0; i < h; i++){
    for (int j = 0; j < w; j++){
      tmp_im.data[i*w+j] = !!im->data[i*2*im->w+j*2];
    }
  }
  
  int ws[thick*2+1];
  for (int i = 0; i < thick*2+1; i++){
    float t = (float)i/(float)(thick*2) * 2 - 1;
    float y = sqrt(1-t*t);
    ws[i] = ceil(y*thick);
  }

  thinning_zs(&tmp_im);
  for (int i = 0; i < h; i++){
    for (int j = 0; j < w; j++){
      if (!tmp_im.data[i*w+j]){
        continue;
      }
      int x0 = j*2;
      int y0 = i*2;
      if (!thick){
        rast_set(im, x0, y0, 0, val);
      }else{
        int yl = MAX(y0-thick,0);
        int yr = MIN(y0+thick+1,im->h);
        for (int i = yl; i < yr; i++){
          int ww = ws[i-y0+thick];
          int xl = MAX(x0-ww,0);
          int xr = MIN(x0+1+ww,im->w);
          memset(im->data + (i * im->w) + xl, val, xr-xl);
        }
      }
    }
  }
}

float polys_closest_points(polyline_t* p0, polyline_t* p1, float* xyxy){
  float md = INFINITY;
  point_t* it = p0->head;
  while (it){
    point_t* jt = p1->head;
    while (jt){
      float d = dist2d(it,jt);
      if (d < md){
        md = d;
        xyxy[0] = it->x;
        xyxy[1] = it->y;
        xyxy[2] = jt->x;
        xyxy[3] = jt->y;
      }
      jt = jt->next;
    }
    it = it->next;
  }
  return md;
}

void draw_bridges(rast_t* im, polyline_t* polys, bbox_t bb, int thick, uint8_t val){
  polyline_t* it = polys;
  while (it){
    if (it->lvl){
      it = it->next;
      continue;
    }
    float line[4];
    float md = INFINITY;
    polyline_t* jt = polys;
    while (jt != it){
      if (jt->lvl){
        jt = jt->next;
        continue;
      }

      float xyxy[4];
      float d = polys_closest_points(it,jt,xyxy);
      if (d < md){
        md = d;
        line[0]=xyxy[0], line[1]=xyxy[1], line[2]=xyxy[2], line[3]=xyxy[3];
      }
      jt = jt->next;
    }
    if (md != INFINITY){
      plot_line(im,w2c_x(line[0]),w2c_y(line[1]),w2c_x(line[2]),w2c_y(line[3]),thick,val);
    }
    it = it->next;
  }
}

void fit_mesh_to_z(mesh_t* mesh, float zh){
  bbox_t bb = bounding_box(mesh);
  float z0 = bb.min.z;
  float z1 = bb.max.z;
  float d = z1-z0;
  float s = zh/d;
  for (int i = 0; i < mesh->n_vtx; i++){
    mesh->vtxs[i].x *= s;
    mesh->vtxs[i].y *= s;
    mesh->vtxs[i].z -= z0;
    mesh->vtxs[i].z *= s;
  }
}

void swap_mesh_yz(mesh_t* mesh){
  for (int i = 0; i < mesh->n_vtx; i++){
    float y = mesh->vtxs[i].y;
    float z = mesh->vtxs[i].z;
    mesh->vtxs[i].x *= -1;
    mesh->vtxs[i].y = z;
    mesh->vtxs[i].z = y;
  }
}


polyline_t* search_iden_seg(polyline_t* segs, polyline_t* pts){
  polyline_t* it = segs;

  while (it){
    if (dist2d(it->head, pts->head) == 0 && dist2d(it->tail, pts->tail) == 0){
      return it;
    }
    if (dist2d(it->head, pts->tail) == 0 && dist2d(it->tail, pts->head) == 0){
      return it;
    }
    it = it->next;
  }
  return NULL;
}


void floodfill(rast_t* im){
  polyline_t* Q = new_polyline();
  add_point_to_polyline(Q,0,0);
  rast_set(im,0,0,0,255);

  while (Q->size > 0){
    point_t* n = Q->head;
    Q->head = n->next;
    Q->size--;

    // printf("%d %f %f\n",Q->size,n->x,n->y);
    
    if (!rast_get_white_border(im,n->x-1,n->y,0)){
      add_point_to_polyline(Q, n->x-1,n->y);
      rast_set(im,n->x-1,n->y,0,255);
    }
    if (!rast_get_white_border(im,n->x+1,n->y,0)){
      add_point_to_polyline(Q, n->x+1,n->y);
      rast_set(im,n->x+1,n->y,0,255);
    }
    if (!rast_get_white_border(im,n->x,n->y-1,0)){
      add_point_to_polyline(Q, n->x,n->y-1);
      rast_set(im,n->x,n->y-1,0,255);
    }
    if (!rast_get_white_border(im,n->x,n->y+1,0)){
      add_point_to_polyline(Q, n->x,n->y+1);
      rast_set(im,n->x,n->y+1,0,255);
    }
    free(n);
  }
  destroy_polylines(Q);
}

rast_t mesh_to_voxels(mesh_t* mesh, float step, int min_thick, int min_thin, int no_hole, int bridge, int do_blur, int check_dup, int inv_flood, const char* dump_pth){
  static char fname[1024];

  bbox_t bb = bounding_box(mesh);
  float z0 = bb.min.z;
  float z1 = bb.max.z;

  float pad = step*8;
  bb.min.x -= pad;
  bb.min.y -= pad;
  bb.max.x += pad;
  bb.max.y += pad;

  float mw = bb.max.x - bb.min.x;
  float mh = bb.max.y - bb.min.y;
  float md = bb.max.z - bb.min.z;

  int d = md / step;
  int w = (md/step) * (mw/md);
  int h = (md/step) * (mh/md);

  rast_t vox = new_rast(w,h,d);

  rast_t tmp0 = new_rast(w,h,1);

  float* zlims = malloc(sizeof(float)*mesh->n_tri*2);
  for (int i = 0; i < mesh->n_tri; i++){
    vec3_t p0 = mesh->vtxs[mesh->tris[i*3]];
    vec3_t p1 = mesh->vtxs[mesh->tris[i*3+1]];
    vec3_t p2 = mesh->vtxs[mesh->tris[i*3+2]];
    zlims[i*2] =   fmin(fmin(p0.z,p1.z),p2.z);
    zlims[i*2+1] = fmax(fmax(p0.z,p1.z),p2.z);
  }

  // z = 100;
  // z = 100.973221;
  for (int zi = 0; zi < d; zi++){
    float z = z0+(float)zi*step+ step*0.01;

    // z = 8.673220;
    // z = 49.373222;
    // z = 53.173222;
    // z = 69.173225;
    // z = 0.398000;
    // z = 86.298004;
    // z = 9.298000;
    // z = 18.798000;
    if (zi % 10 == 0) printf("[voxelize] working on z=%f (%d/%d)\n",z,zi,d);
    polyline_t* segs = NULL;

    for (int i = 0; i < mesh->n_tri; i++){
      // printf("%d/%d\n",i,mesh->n_tri);
      if (z < zlims[i*2] || z >= zlims[i*2+1]){
        continue;
      }

      vec3_t p0 = mesh->vtxs[mesh->tris[i*3+0]];
      vec3_t p1 = mesh->vtxs[mesh->tris[i*3+1]];
      vec3_t p2 = mesh->vtxs[mesh->tris[i*3+2]];
      vec3_t pt;
      polyline_t* pts = new_polyline();

      if (seg_isect_z_plane(z,p0,p1,&pt)){
        add_point_to_polyline(pts, pt.x,pt.y);
      }
      if (seg_isect_z_plane(z,p1,p2,&pt)){
        add_point_to_polyline(pts, pt.x,pt.y);
      }
      if (seg_isect_z_plane(z,p2,p0,&pt)){
        add_point_to_polyline(pts, pt.x,pt.y);
      }
      if (dist2d(pts->head,pts->tail)>0){
        if (!check_dup || !search_iden_seg(segs,pts)){
          segs = prepend_polyline(segs,pts);
        }else{
          destroy_polylines(pts);
        }
      }else{
        // printf("%f,%f,%f %f,%f,%f %f,%f,%f %f,%f\n",p0.x,p0.y,p0.z,p1.x,p1.y,p1.z,p2.x,p2.y,p2.z,pt.x,pt.y);
        destroy_polylines(pts);
      }
    }

    if (dump_pth){
      sprintf(fname,"%s/%.3f_segm.svg",dump_pth,z);
      polylines_to_svg(fname,bb.min.x,bb.min.y,bb.max.x,bb.max.y,segs);
    }

    rast_t im;
    im.w = w;
    im.h = h;
    im.d = 1;
    im.data = vox.data + ((uint64_t)zi*(uint64_t)w*(uint64_t)h);

    if (inv_flood){
      int lw = inv_flood-1;

      draw_polygon_outlines(&im,segs,bb,-1,lw,255);

      if (dump_pth){
        sprintf(fname,"%s/%.03f_outl.gif",dump_pth,z);
        write_gif(fname,&im);
      }

      floodfill(&im);
      draw_polygon_outlines(&im,segs,bb,-1,lw,0);

      for (int i = 0; i < w*h; i++){
        im.data[i] = 255 - im.data[i];
      }

    }else{

      polyline_t* polys = connect_segs(segs,0);
      // polyline_t* polys = segs;

      if (dump_pth){
        sprintf(fname,"%s/%.3f_poly.svg",dump_pth,z);
        polylines_to_svg(fname,bb.min.x,bb.min.y,bb.max.x,bb.max.y,polys);
      }

      int nlvl = polys_check_internal(polys)+1;

      if (!min_thin && !min_thick){
        draw_polygon_fills(&tmp0,polys,bb,no_hole?0:-1,255);

      }else{
        for (int k = 0; k < nlvl; k++){
          if (k % 2 == 0){
            rast_clear(&tmp0);
            draw_polygon_fills(&tmp0,polys,bb,k,255);
            rast_blend(&im, &tmp0, BLEND_MAX);
          }else{
            draw_polygon_fills(&im,polys,bb,k,0);

            if (min_thick){
              draw_polygon_outlines(&im,polys,bb,k-1,min_thick,255);
              rast_blend(&im, &tmp0, BLEND_MIN);
            }
          }
          if (no_hole) break;
        }
      }

      if (bridge){
        draw_bridges(&im,polys,bb,bridge,255);
      }

      destroy_polylines(polys);
    }

    if (min_thin){
      draw_skeleton(&im,min_thin,255);
    }

    if (do_blur){
      rast_clear(&tmp0);
      for (int i = 1; i < im.h-1; i++){
        for (int j = 1; j < im.w-1; j++){
          int a = rast_get(&im,j-1,i,0);
          int b = rast_get(&im,j  ,i,0);
          int c = rast_get(&im,j+1,i,0);
          rast_set(&tmp0,j,i,0,a*0.25+b*0.5+c*0.25);
        }
      }
      for (int i = 1; i < im.h-1; i++){
        for (int j = 1; j < im.w-1; j++){
          int a = rast_get(&tmp0,j,i-1,0);
          int b = rast_get(&tmp0,j,i  ,0);
          int c = rast_get(&tmp0,j,i+1,0);
          rast_set(&im,j,i,0,a*0.25+b*0.5+c*0.25);
        }
      }
    }

    if (dump_pth){
      sprintf(fname,"%s/%.03f_rast.gif",dump_pth,z);
      write_gif(fname,&im);
    }

    
    // break;
  }

  rast_destroy(&tmp0);
  free(zlims);
  return vox;

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
  // for (int i = 0; i < n; i++){
  //   fputc(vox->data[i],fp);
  // }
  fclose(fp);
}

void print_help(){
  printf("turn mesh to voxels, making printability edits at raster level              \n");
  printf("USAGE: mesh_to_voxels [options] mesh.obj                                    \n");
  printf("OPTIONS:                                                                    \n");
  printf("--output, -o   path/out.bin    output path, voxel binary format: (LE)       \n");
  printf("                               w (u32), h (u32), d (u32), data (u8,u8,u8,...\n");
  printf("                               voxel(x,y,z) = data[z*w*h + y*w + x]         \n");
  printf("--step         0.1             slicing z step (layer height)                \n");
  printf("--min_thick    8               minimum shell thickness (px)                 \n");
  printf("                               (remove shell-scratching internal holes)     \n");
  printf("--min_thin     0               minimum feature thinness (px)                \n");
  printf("                               (increase size of tiny/skinny parts)         \n");
  printf("--no_hole                      flag: get rid of all holes                   \n");
  printf("--bridge       0               bridge all shells given bridge thickness     \n");
  printf("                               (for better vase, also removes holes,        \n");
  printf("                               default 0 = no bridging)                     \n");
  printf("--swap_yz                      swap y and z axes of mesh                    \n");
  printf("--fit_z        100             pre-scale mesh to given z dimension          \n");
  printf("                               (default 0 = no scaling)                     \n");
  printf("--blur                         flag: add slight blur (simulate anti-alias,  \n");
  printf("                               xy only, might improve marching cubes)       \n");
  printf("--check_dup                    flag: check duplicate segments               \n");
  printf("--inv_flood    0               use inverted flood-fill to generate contours \n");
  printf("                               (more robust, invalidate most other features,\n");
  printf("                               for truly nefarious meshes only,             \n");
  printf("                               number = dilation amount, 0 = disabled)      \n");
  printf("--dump         path/folder     dump intermediate visuals (for debugging)    \n");
  printf("--help                         print this message                           \n");
}

int main(int argc, char** argv){

  const char* dump_pth = NULL;
  const char* outp_pth = NULL;
  const char* inpt_pth = NULL;
  int min_thick = 8;
  int min_thin = 0;
  float step = 0.1;
  float fit_z = 0;
  int swap_yz = 0;
  int no_hole = 0;
  int bridge = 0;
  int do_blur = 0;
  int check_dup = 0;
  int inv_flood = 0;
  
  int i = 1;
  while (i < argc){
    if (!strcmp(argv[i],"--min_thick")){
      min_thick = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--min_thin")){
      min_thin = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--bridge")){
      bridge = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--no_hole")){
      no_hole = 1;
      i+=1;
    }else if (!strcmp(argv[i],"--dump")){
      dump_pth = argv[i+1];
      i+=2;
    }else if (!strcmp(argv[i],"--output") || !strcmp(argv[i],"-o")){
      outp_pth = argv[i+1];
      i+=2;
    }else if (!strcmp(argv[i],"--fit_z")){
      fit_z = atof(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--swap_yz")){
      swap_yz = 1;
      i++;
    }else if (!strcmp(argv[i],"--step")){
      step = atof(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--blur")){
      do_blur = 1;
      i++;
    }else if (!strcmp(argv[i],"--check_dup")){
      check_dup = 1;
      i++;
    }else if (!strcmp(argv[i],"--inv_flood")){
      inv_flood = atoi(argv[i+1]);
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
  
  mesh_t mesh = load_obj(inpt_pth);
  printf("[read obj] model loaded.\n");
  
  if (swap_yz){
    swap_mesh_yz(&mesh);  
  }

  if (fit_z > 0){
    fit_mesh_to_z(&mesh,fit_z);  
  }

  rast_t vox = mesh_to_voxels(&mesh,step,min_thick,min_thin,no_hole,bridge,do_blur,check_dup,inv_flood,dump_pth);
  printf("[write bin] generating output...\n");
  write_voxels_to_file(outp_pth, &vox);

  rast_destroy(&vox);
  free(mesh.vtxs);
  free(mesh.tris);
  return 0;
}