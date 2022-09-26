#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#define SOLI 255
#define SPRT 2
#define VOID 0

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
uint8_t rast_set(rast_t* r, int x, int y, int z, uint8_t v){
  if (x < 0 || x >= r->w || y < 0 || y >= r->h || z < 0 || z >= r->d){
    return 0;
  }
  return r->data[z*(r->w*r->h)+y*(r->w)+x] = v;
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


void erode(rast_t* r,rast_t* r1, int k){
  for (int z = 0; z < r->d; z++){
    for (int i = 0; i < r->h; i++){
      for (int j = 0; j < r->w; j++){
        int v = rast_get(r,j,i,z);
        rast_set(r1,j,i,z,v);
        if (v == SPRT){
          for (int ki = -k; ki <= k; ki++){
            if (rast_get(r,j+ki,i,z)!=SPRT){
              rast_set(r1,j,i,z,VOID);
              goto nextpix0;
            }
          }
        }
        nextpix0:
        continue;
      }
    }
  }
  for (int z = 0; z < r->d; z++){
    for (int i = 0; i < r->h; i++){
      for (int j = 0; j < r->w; j++){
        int v = rast_get(r1,j,i,z);
        rast_set(r,j,i,z,v);
        if (v == SPRT){
          for (int ki = -k; ki <= k; ki++){
            if (rast_get(r1,j,i+ki,z)!=SPRT){
              rast_set(r,j,i,z,VOID);
              goto nextpix1;
            }
          }
        }
        nextpix1:
        continue;
      }
    }
  }
}


void dilate(rast_t* r,rast_t* r1, int k, int kz){
  printf("[dilate] x...\n");
  for (int z = 0; z < r->d; z++){
    for (int i = 0; i < r->h; i++){
      for (int j = 0; j < r->w; j++){
        int v = rast_get(r,j,i,z);
        rast_set(r1,j,i,z,v);
        if (v == VOID){
          for (int ki = -k; ki <= k; ki++){
            if (rast_get(r,j+ki,i,z)==SPRT){
              rast_set(r1,j,i,z,SPRT);
              goto nextpix0;
            }
          }
        }
        nextpix0:
        continue;
      }
    }
  }
  printf("[dilate] y...\n");
  for (int z = 0; z < r->d; z++){
    for (int i = 0; i < r->h; i++){
      for (int j = 0; j < r->w; j++){
        int v = rast_get(r1,j,i,z);
        rast_set(r,j,i,z,v);
        if (v == VOID){
          for (int ki = -k; ki <= k; ki++){
            if (rast_get(r1,j,i+ki,z)==SPRT){
              rast_set(r,j,i,z,SPRT);
              goto nextpix1;
            }
          }
        }
        nextpix1:
        continue;
      }
    }
  }
  printf("[dilate] z...\n");
  for (int z = 0; z < r->d; z++){
    for (int i = 0; i < r->h; i++){
      for (int j = 0; j < r->w; j++){
        int v = rast_get(r,j,i,z);
        rast_set(r1,j,i,z,v);
        if (v == VOID){
          for (int ki = -kz; ki <= kz; ki++){
            if (rast_get(r,j,i,z+ki)==SPRT){
              rast_set(r1,j,i,z,SPRT);
              goto nextpix2;
            }
          }
        }
        nextpix2:
        continue;
      }
    }
  }
  memcpy(r->data,r1->data,r->w*r->h*r->d*sizeof(uint8_t));
}

int hang_cot = 4;
int ero_xy = 1;
int dil_xy = 24;
int dil_z = 16;
int sink = 0;

void generate_support(rast_t* vox){
  rast_t tmp0 = new_rast(vox->w,vox->h,vox->d);
  for (int i = 0; i < vox->w*vox->h*vox->d; i++){
    vox->data[i] = vox->data[i] > 128 ? SOLI : VOID;
  }
  printf("supporting...\n");
  for (int z = vox->d-1; z >= 0; z--){
    for (int y = 0; y < vox->h; y++){
      for (int x = 0; x < vox->w; x++){
        if (!rast_get(vox,x,y,z)){
          
          int p = rast_get(vox,x,y,z+1);
          if (p == SOLI){
            int hang = 1;
            for (int i = -hang_cot; i <= hang_cot; i++){
              for (int j = -hang_cot; j <= hang_cot; j++){
                if (rast_get(vox,x+j,y+i,z-1)){
                  hang = 0;
                  goto done_check_hang;
                }
              }
            }
            done_check_hang:
            if (hang){
              rast_set(vox,x,y,z,SPRT);
            }
          }else if (p == SPRT){
            rast_set(vox,x,y,z,SPRT);
          }
        }
      }
    }
  }
  if (ero_xy){
    printf("eroding...\n");
    erode(vox,&tmp0,ero_xy);
  }
  printf("dilating...\n");
  dilate(vox,&tmp0,dil_xy,dil_z);

  printf("supporting support...\n");
  for (int z = vox->d-1; z >= 0; z--){
    for (int y = 0; y < vox->h; y++){
      for (int x = 0; x < vox->w; x++){
        if (!rast_get(vox,x,y,z)){          
          int p = rast_get(vox,x,y,z+1);
          if (p == SPRT){
            rast_set(vox,x,y,z,SPRT);
          }
        }
      }
    }
  }

  rast_destroy(&tmp0);
}


void print_help(){
  printf("generate 3d printing support for voxels. 0=void, 255=solid, 2=support.      \n");
  printf("USAGE: generate_support [options] input.bin                                 \n");
  printf("OPTIONS:                                                                    \n");
  printf("--output, -o   path/out.bin    output path, voxel binary format: (BE)       \n");
  printf("                               w (u32), h (u32), d (u32), data (u8,u8,u8,...\n");
  printf("                               voxel(x,y,z) = data[z*w*h + y*w + x]         \n");
  printf("--hang         4               max overhang slope cotangent                 \n");
  printf("                               (larger the number, less the support)        \n");
  printf("--ero_xy       1               erosion on xy plane (remove tiny supports)   \n");
  printf("--dil_xy       24              dilation on xy plane (merge nearby supports) \n");
  printf("--dil_z        16              dilation on z axis (envelop bottom of hangs) \n"); 
  printf("--sink         0               number of layers to sink into build plane (z)\n");  
}

int main(int argc, char** argv){

  const char* outp_pth = NULL;
  const char* inpt_pth = NULL;

  int i = 1;
  while (i < argc){
    if (!strcmp(argv[i],"--hang")){
      hang_cot = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--ero_xy")){
      ero_xy = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--dil_xy")){
      dil_xy = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--dil_z")){
      dil_z = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--sink")){
      sink = atoi(argv[i+1]);
      i+=2;
    }else if (!strcmp(argv[i],"--output") || !strcmp(argv[i],"-o")){
      outp_pth = argv[i+1];
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

  rast_t vox = read_voxels_from_file(inpt_pth);
  printf("[voxel] loaded.\n");

  rast_t vox2;
  vox2.w = vox.w;
  vox2.h = vox.h;
  vox2.d = vox.d-sink;
  vox2.data = vox.data + ((uint64_t)sink*(uint64_t)vox.w*(uint64_t)vox.h);

  generate_support(&vox2);
  printf("[voxel] writing...\n");
  write_voxels_to_file(outp_pth,&vox2);
  return 0;
}