// https://en.wikipedia.org/wiki/Boy%27s_surface

#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>

typedef struct {
  uint8_t* data;
  int w;
  int h;
  int d;
  int cap;
} rast_t;

rast_t new_rast(int w, int h, int d){
  rast_t r;
  r.cap = w*h*d;
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

float dist(float x0, float y0, float z0, float x1, float y1, float z1){
  float dx = x1-x0;
  float dy = y1-y0;
  float dz = z1-z0;
  return sqrt(dx*dx+dy*dy+dz*dz);
}

void setball(rast_t* vox, int x, int y, int z, int r){
  for (int i = z-r; i <= z+r; i++){
    for (int j = y-r; j<= y+r; j++){
      for (int k = x-r; k<= x+r; k++){
        if (dist(k,j,i,x,y,z) <= r+1){
          rast_set(vox,k,j,i,255);
        }
      }
    }
  }
}

const int mQ = 100000000;
int Q[mQ];
int zQ = 0;
int nQ = 0;

int Q_len(){
  if (zQ > nQ){
    return nQ+mQ-zQ;
  }
  return nQ-zQ;
}
int Q_pop(){
  int v = Q[zQ];
  zQ = (zQ+1)%mQ;
  return v;
}
void Q_push(int v){
  Q[nQ] = v;
  nQ = (nQ+1)%mQ;
}

void floodfill3d(rast_t* vox){
  nQ = 0;
  zQ = 0;
  Q_push(0);
  rast_set(vox,0,0,0,1);

  while (Q_len()){
    int p = Q_pop();

    int z = p / (vox->w*vox->h);
    int y = (p - (z*vox->w*vox->h)) / vox->w;
    int x = p % vox->w;

    if (!rast_get_white_border(vox,x-1,y,z)){
      Q_push(p-1);
      rast_set(vox,x-1,y,z,1);
    }
    if (!rast_get_white_border(vox,x+1,y,z)){
      Q_push(p+1);
      rast_set(vox,x+1,y,z,1);
    }
    if (!rast_get_white_border(vox,x,y-1,z)){
      Q_push(p-vox->w);
      rast_set(vox,x,y-1,z,1);
    }
    if (!rast_get_white_border(vox,x,y+1,z)){
      Q_push(p+vox->w);
      rast_set(vox,x,y+1,z,1);
    }
    if (!rast_get_white_border(vox,x,y,z-1)){
      Q_push(p-vox->w*vox->h);
      rast_set(vox,x,y,z-1,1);
    }
    if (!rast_get_white_border(vox,x,y,z+1)){
      Q_push(p+vox->w*vox->h);
      rast_set(vox,x,y,z+1,1);
    }
  }
  printf("%d %d\n",zQ,nQ);
}


int solid = 1;

int main(){
  int n = 2000;
  int r = 5;

  rast_t vox = new_rast(600+r*2,600+r*2,600+r*2);

  float xmin = INFINITY;
  float xmax = -INFINITY;
  float ymin = INFINITY;
  float ymax = -INFINITY;
  float zmin = INFINITY;
  float zmax = -INFINITY;


  xmin = -1.333879;
  xmax = 1.265682;
  ymin = -1.333483;
  ymax = 1.255926;
  zmin = -0.442630;
  zmax = 2.000000;

  for (int i = 0; i < n; i++){
    if (i%100==0)printf("%d/%d\n",i,n);
    for (int j = 0; j < n; j++){
      double complex vv = (((float)i/(float)n-0.5)*2) * I + (((float)j/(float)n-0.5)*2) ;
      // printf("%lf %lf\n",creal(vv),cimag(vv));
      if (cabs(vv)>1){
        continue;
      }
      double complex vv3 = vv*vv*vv;
      double complex vv4 = vv3*vv;
      double complex vv6 = vv4*vv*vv;
      double g1 = -1.5 * cimag(( vv * (1 - vv4) ) / (vv6 + sqrt(5) * vv3 - 1));
      double g2 = -1.5 * creal(( vv * (1 + vv4) ) / (vv6 + sqrt(5) * vv3 - 1));
      double g3 = -0.5 + cimag((1 + vv6)  / (vv6 + sqrt(5) * vv3 - 1));
      double gg =  1.0/(g1*g1+g2*g2+g3*g3);
      g1*=gg;
      g2*=gg;
      g3*=-gg;

      // xmin = fmin(xmin,g1); xmax = fmax(xmax,g1);
      // ymin = fmin(ymin,g2); ymax = fmax(ymax,g2);
      // zmin = fmin(zmin,g3); zmax = fmax(zmax,g3);

      int x = ( (g1 - xmin)/(xmax-xmin) ) * (600-r*2) + r;
      int y = ( (g2 - ymin)/(ymax-ymin) ) * (600-r*2) + r;
      int z = ( (g3 - zmin)/(zmax-zmin) ) * (600-r*2) + r;

      
      setball(&vox,x,y,z,r);
    }
  }

  if (solid){
    floodfill3d(&vox);
    for (int i = 0; i < vox.w*vox.h*vox.d; i++){
      vox.data[i] = (vox.data[i] == 0 || vox.data[i] == 255) ? 255 : 0;
    }
  }    

  // printf("%f %f %f %f %f %f\n",xmin,xmax, ymin,ymax, zmin,zmax);

  char fname[128];
  sprintf(fname,"output/boy-%s.bin",solid?"solid":"surf");
  write_voxels_to_file(fname,&vox);
}

