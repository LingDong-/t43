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

void export_numpy(const char* pth, rast_t* r){
  FILE *fp;
  fp = fopen(pth, "w+");
  fputc(0x93,fp);
  fputs("NUMPY",fp);
  fputc(0x1,fp);
  fputc(0x0,fp);
  fputc(118,fp);
  fputc(0,fp);
  fprintf(fp,"{'descr':'|u1','fortran_order':False,'shape':(%9d,%9d,%9d)}\n",r->d,r->h,r->w);
  fputs("                                        ",fp); //40 spaces
  int n = r->w*r->h*r->d;
  fwrite((r->data),1,n,fp);
  fclose(fp);
}

#define NPY_READ_INFO_TIL(c,yesno)  while((info[i]==(c))^yesno){i++;if(i>=hlen||info[i]=='\n'){printf("[error] file read failure: cannot find expected symbol.\n");exit(1);}}

rast_t import_numpy(const char* pth){
  FILE * fp;
  fp = fopen(pth,"r");
  if (!fp){ printf("cannot open file."); exit(1); }
  char magic[9];
  fread(&magic, 1, 8, fp);
  magic[8] = 0;
  if (strcmp(magic,"\x93NUMPY\1\0")){
    printf("[warn] magic number/version mismatch. (only numpy 1.0 supported).\n");
  }
  int hlen = fgetc(fp);
  hlen = (fgetc(fp)<<8) | hlen;

  char info[hlen];
  fread(&info, 1, hlen, fp);
  int i = 4;
  int w, h, d;
  while (1){
    if (info[i-4] == 'd' && info[i-3] == 'e' && info[i-2] == 's' && info[i-1] == 'c' && info[i] == 'r'){
      NPY_READ_INFO_TIL(':',1); i++;
      NPY_READ_INFO_TIL(' ',0);
      i++;
      if (info[i] != '|' || info[i+1] != 'u' || info[i+2] != '1'){
        printf("[error] only dtype='|u1' (unsigned byte) supported! (found '%c%c%c')\n",info[i],info[i+1],info[i+2]);
        exit(1);
      }
    }else if (info[i-4] == 'o' && info[i-3] == 'r' && info[i-2] == 'd' && info[i-1] == 'e' && info[i] == 'r'){
      NPY_READ_INFO_TIL(':',1); i++;
      NPY_READ_INFO_TIL(' ',0);
      if (info[i] != 'F' && info[i] != '0'){
        printf("[warn] fortran_order seems to be True, but only C order is supported!\n");
      }
    }else if (info[i-4] == 's' && info[i-3] == 'h' && info[i-2] == 'a' && info[i-1] == 'p' && info[i] == 'e'){
      NPY_READ_INFO_TIL('(',1); i++;
      d = atoi(info+i);
      NPY_READ_INFO_TIL(',',1); i++;
      h = atoi(info+i);
      NPY_READ_INFO_TIL(',',1); i++;
      w = atoi(info+i);
    }
    i++;
    if (i > hlen || info[i] == '\n'){
      break;
    }
  }
  int n = w*h*d;
  rast_t r = new_rast(w,h,d);
  fread(r.data,1,n,fp);
  return r;
}


void export_runlen(const char* pth, rast_t* r){
  FILE *fp;
  if (pth){
    fp = fopen(pth, "wb");
  }else{
    fp = stdout;
  }

  uint16_t n = 0;
  int c = 0;
  fwrite(&(r->w),4,1,fp);
  fwrite(&(r->h),4,1,fp);
  fwrite(&(r->d),4,1,fp);
  for (int i = 0; i < r->w*r->h*r->d; i++){
    if (!!(r->data[i]) != c){
      fwrite(&n,2,1,fp);
      c ^= 1;
      n = 1;
    }else{
      n++;
      if (n >= 65535){
        fwrite(&n,2,1,fp);
        n = 0;
        fwrite(&n,2,1,fp);
      }
    }
  }
  if (n){
    fwrite(&n,2,1,fp);
  }
  fclose(fp);
}

rast_t import_runlen(const char* pth){
  FILE * fp;
  fp = fopen(pth,"r");
  if (!fp){ printf("cannot open file."); exit(1); }
  int w, h, d;
  fread(&w, 4, 1, fp);
  fread(&h, 4, 1, fp);
  fread(&d, 4, 1, fp);
  int n = w*h*d;
  rast_t r = new_rast(w,h,d);
  int l;
  int c = 0;
  int i = 0;
  while (fread(&l,2,1,fp)){
    for (int j = 0; j < l; j++){
      r.data[i++] = c;
    }
    c = 255-c;
  }
  return r;
}


void print_help(){
  printf("convert voxel binary format used by this toolkit from/to common formats     \n");
  printf("currently supported formats:                                                \n");
  printf("*.bin: our format:   w (u32 LE), h (u32), d (u32), data (u8,u8,u8,...       \n");
  printf("                     voxel(x,y,z) = data[z*w*h + y*w + x]                   \n");
  printf("*.npy: numpy         only supports dtype='|u1' (unsigned bytes), C-order    \n");
  printf("                     shape=(d,h,w). see:                                    \n");
  printf("     https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html \n");
  printf("                     online npy visualizer https://tensor-viewer.glitch.me/ \n");
  printf("*.rle.bin: run-      w (u32 LE), h (u32), d (u32),                          \n");
  printf(" length encoding     # 0's (u16),  # 255's (u16), # 0's, # 255's, # 0's, ...\n");
  printf("                     (more space-efficient for storage)                     \n");
  printf("USAGE: voxel_convert [options] input.ext                                    \n");
  printf("OPTIONS:                                                                    \n");
  printf("--output, -o   path/out.ext    output path                                  \n");
  printf("--from         npy             npy, bin, rle                                \n");
  printf("--to           bin             npy, bin, rle                                \n");
  printf("--help                         print this message                           \n");
}

int main(int argc, char** argv){

  const char* outp_pth = NULL;
  const char* inpt_pth = NULL;
  const char* from = "npy";
  const char* to = "bin";

  int i = 1;
  while (i < argc){
    if (!strcmp(argv[i],"--output") || !strcmp(argv[i],"-o")){
      outp_pth = argv[i+1];
      i+=2;
    }else if (!strcmp(argv[i],"--from")){
      from = argv[i+1];
      i+=2;
    }else if (!strcmp(argv[i],"--to")){
      to = argv[i+1];
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

  rast_t vox;
  if (!strcmp(from,"bin")){
    vox = read_voxels_from_file(inpt_pth);
  }else if (!strcmp(from,"npy")){
    vox = import_numpy(inpt_pth);
  }else if (!strcmp(from,"rle")){
    vox = import_runlen(inpt_pth);
  }
  printf("%s read: %d x %d x %d\n",from,vox.w,vox.h,vox.d);

  if (!strcmp(to,"bin")){
    write_voxels_to_file(outp_pth,&vox);
  }else if (!strcmp(to,"npy")){
    export_numpy(outp_pth,&vox);
  }else if (!strcmp(to,"rle")){
    export_runlen(outp_pth,&vox);
  }

  printf("%s written.\n",to);
}