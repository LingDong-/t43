// https://en.wikipedia.org/wiki/Trefoil_knot

#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

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

int p = 3;
int q = 5;

int size = 120;

int main(int argc, char** argv){
  if (argc > 1){
    p = atoi(argv[1]);
  }
  if (argc > 2){
    q = atoi(argv[2]);
  }
  int n = 1500;
  int w = 60;

  rast_t vox = new_rast(size*6+w*2,size*6+w*2,size*2+w*2);

  for (int i = 0; i < n; i++){
    if (i%100==0)printf("%d/%d\n",i,n);
    float t = (float)i/n * M_PI * 2;

    int x = (int)( (size*3+w)+size*((2+cos(q*t))*cos(p*t)) );
    int y = (int)( (size*3+w)+size*((2+cos(q*t))*sin(p*t)) );
    int z = (int)( (size+w)+size*(sin(q*t) ) );
    
    setball(&vox,x,y,z,w);
  }
  char fname[128];
  sprintf(fname,"output/knot%d%d.bin",p,q);
  write_voxels_to_file(fname,&vox);


}