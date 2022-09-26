#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

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
uint8_t rast_get_cp_border(rast_t* r, int x, int y, int z){
  x = MIN(MAX(x,0),r->w-1);
  y = MIN(MAX(y,0),r->h-1);
  z = MIN(MAX(z,0),r->d-1);
  return r->data[z*(r->w*r->h)+y*(r->w)+x];
}
uint8_t rast_set(rast_t* r, int x, int y, int z, uint8_t v){
  if (x < 0 || x >= r->w || y < 0 || y >= r->h || z < 0 || z >= r->d){
    return 0;
  }
  return r->data[z*(r->w*r->h)+y*(r->w)+x] = v;
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

void blur_voxels(rast_t* vox, int rad, int no_z, const char* dump_pth){
  static char fname[1024];

  int ksize = rad*2+1;
  float kern[ksize];

  rast_t voy = new_rast(vox->w,vox->h,vox->d);

  float sigma = 0.3*((float)(ksize-1)*0.5 - 1) + 0.8;
  float ss2 = sigma*sigma*2;
  int i; for (i = 0; i < ksize; i++){
    float x = (float)(i-ksize/2);
    float z = exp(-(x*x)/(ss2))/(2.5066282746*sigma);
    kern[i] = z;
  }

  for (int z = 0; z < vox->d; z++){
    if (z % 100 == 0) printf("[blur] x axis, z = %d / %d\n",z,vox->d);
    for (int y = 0; y < vox->h; y++ ){
      for (int x= 0; x < vox->w; x++ ){
        float sum = 0;
        for (int k = 0; k < ksize; k++ ){
          sum += rast_get_cp_border(vox,x-rad+k,y,z) * kern[k];
        }
        rast_set(&voy,x,y,z,sum);
      }
    }
  }
  for (int z = 0; z < vox->d; z++){
    if (z % 100 == 0) printf("[blur] y axis, z = %d / %d\n",z,vox->d);
    for (int y = 0; y < vox->h; y++ ){
      for (int x= 0; x < vox->w; x++ ){
        float sum = 0;
        for (int k = 0; k < ksize; k++ ){
          sum += rast_get_cp_border(&voy,x,y-rad+k,z) * kern[k];
        }
        rast_set(vox,x,y,z,sum);
      }
    }
    if (no_z && dump_pth){
      rast_t im;
      im.w = voy.w;
      im.h = voy.h;
      im.d = 1;
      im.data = voy.data + ((uint64_t)z*(uint64_t)voy.w*(uint64_t)voy.h);
      sprintf(fname,"%s/%03d_blur.gif",dump_pth,z);
      write_gif(fname,&im);
    }
  }

  if (!no_z){
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
  }

  memcpy(vox->data,voy.data,vox->w*vox->h*vox->d);
  free(voy.data);

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
  // for (int i = 0; i < n; i++){
  //   fputc(vox->data[i],fp);
  // }
  fclose(fp);
}


void print_help(){
  printf("gaussian blur voxels, for cube-marching smoothness                          \n");
  printf("USAGE: blur_voxels [options] input.bin                                      \n");
  printf("OPTIONS:                                                                    \n");
  printf("--output, -o   path/out.bin    output path, voxel binary format: (BE)       \n");
  printf("                               w (u32), h (u32), d (u32), data (u8,u8,u8,...\n");
  printf("                               voxel(x,y,z) = data[z*w*h + y*w + x]         \n");
  printf("--radius, -r   2               blur radius (x2+1=kernel size)               \n");
  printf("--no_z                         don't blur on z axis                         \n");
  printf("--dump         path/folder     dump intermediate visuals (for debugging)    \n");
  printf("--help                         print this message                           \n");
}

  
int main(int argc, const char** argv){
  const char* dump_pth = NULL;
  const char* outp_pth = NULL;
  const char* inpt_pth = NULL;
  int rad = 2;
  int no_z = 0;

  int i = 1;
  while (i < argc){
    if (!strcmp(argv[i],"--output") || !strcmp(argv[i],"-o")){
      outp_pth = argv[i+1];
      i+=2;
    }else if (!strcmp(argv[i],"--radius") || !strcmp(argv[i],"-r")){
      rad = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--no_z")){
      no_z = 1;
      i++;
    }else if (!strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")){
      print_help();
      exit(0);
    }else if (!strcmp(argv[i],"--dump")){
      dump_pth = argv[i+1];
      i+=2;
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

  rast_t vox = read_voxels_from_file(inpt_pth);
  printf("[read bin] voxel loaded.\n");
  blur_voxels(&vox,rad,no_z,dump_pth);

  printf("[write bin] generating output...\n");
  write_voxels_to_file(outp_pth, &vox);
  free(vox.data);

  return 0;
}