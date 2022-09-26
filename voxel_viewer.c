#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define V3_CROSS(a1,a2,a3,b1,b2,b3) {(a2)*(b3)-(a3)*(b2),(a3)*(b1)-(a1)*(b3),(a1)*(b2)-(a2)*(b1)}
#define V3_DOT(a1,a2,a3,b1,b2,b3)   ((a1)*(b1)+(a2)*(b2)+(a3)*(b3))
#define V3_MAG(a1,a2,a3) (sqrt((a1)*(a1)+(a2)*(a2)+(a3)*(a3)))

#define MAT_IDEN {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1}
#define MAT_ROTX(a) {1,0,0,0, 0,cos(a),sin(a),0, 0,-sin(a),cos(a),0, 0,0,0,1}
#define MAT_ROTY(a) {cos(a),0,-sin(a),0, 0,1,0,0, sin(a),0,cos(a),0, 0,0,0,1}
#define MAT_ROTZ(a) {cos(a),sin(a),0,0,-sin(a),cos(a),0,0, 0,0,1,0, 0,0,0,1}
#define MAT_TRSL(x,y,z) {1,0,0,0, 0,1,0,0, 0,0,1,0, x,y,z,1}
#define MAT_SCAL(x,y,z) {x,0,0,0, 0,y,0,0, 0,0,z,0, 0,0,0,1}
#define MAT_MULT(A,B) {(A)[0]*(B)[0]+(A)[1]*(B)[4]+(A)[2]*(B)[8]+(A)[3]*(B)[12],(A)[0]*(B)[1]+(A)[1]*(B)[5]+(A)[2]*(B)[9]+(A)[3]*(B)[13],(A)[0]*(B)[2]+(A)[1]*(B)[6]+(A)[2]*(B)[10]+(A)[3]*(B)[14],(A)[0]*(B)[3]+(A)[1]*(B)[7]+(A)[2]*(B)[11]+(A)[3]*(B)[15],(A)[4]*(B)[0]+(A)[5]*(B)[4]+(A)[6]*(B)[8]+(A)[7]*(B)[12],(A)[4]*(B)[1]+(A)[5]*(B)[5]+(A)[6]*(B)[9]+(A)[7]*(B)[13],(A)[4]*(B)[2]+(A)[5]*(B)[6]+(A)[6]*(B)[10]+(A)[7]*(B)[14],(A)[4]*(B)[3]+(A)[5]*(B)[7]+(A)[6]*(B)[11]+(A)[7]*(B)[15],(A)[8]*(B)[0]+(A)[9]*(B)[4]+(A)[10]*(B)[8]+(A)[11]*(B)[12],(A)[8]*(B)[1]+(A)[9]*(B)[5]+(A)[10]*(B)[9]+(A)[11]*(B)[13],(A)[8]*(B)[2]+(A)[9]*(B)[6]+(A)[10]*(B)[10]+(A)[11]*(B)[14],(A)[8]*(B)[3]+(A)[9]*(B)[7]+(A)[10]*(B)[11]+(A)[11]*(B)[15],(A)[12]*(B)[0]+(A)[13]*(B)[4]+(A)[14]*(B)[8]+(A)[15]*(B)[12],(A)[12]*(B)[1]+(A)[13]*(B)[5]+(A)[14]*(B)[9]+(A)[15]*(B)[13],(A)[12]*(B)[2]+(A)[13]*(B)[6]+(A)[14]*(B)[10]+(A)[15]*(B)[14],(A)[12]*(B)[3]+(A)[13]*(B)[7]+(A)[14]*(B)[11]+(A)[15]*(B)[15]}

int mode = 0;
int slice = 0;
int curr_slice = 0;
int pix_scale = 1;
int mouse_x = 0;
int mouse_y = 0;
int mouse_is_down = 0;

int colortable [] = {
  0xfff,
  0x0f0,
};

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
  if (x < 0 || x >= r->w || y < 0 || y >= r->h || z < 0 || z >= r->d || z >= slice){
    return 0;
  }
  return r->data[z*(r->w*r->h)+y*(r->w)+x];
}
uint8_t rast_set(rast_t* r, int x, int y, int z, uint8_t v){
  if (x < 0 || x >= r->w || y < 0 || y >= r->h || z < 0 || z >= r->d || z >= slice){
    return 0;
  }
  return r->data[z*(r->w*r->h)+y*(r->w)+x] = v;
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


int W = 1024;
int H = 768;

float rotx = 0;
float roty = 0;
float rotz = 0;
float panx = 0;
float pany = 0;
float panz = -100;
int frame = 0;



rast_t vox;


GLuint dsp;

int shader = 0;

void draw_face(int x,int y,int z,int v,int ax,int ay,int az,int bx,int by,int bz,int cx,int cy,int cz,int dx,int dy,int dz){
  int rgb = colortable[(v+1) % (sizeof(colortable)/sizeof(int))];
  unsigned char r = ((rgb >> 8)  & 0xf)*16;
  unsigned char g = ((rgb >> 4)  & 0xf)*16;
  unsigned char b = ((rgb >> 0)  & 0xf)*16;

  if (shader == 0){
    if (ax == bx && bx == cx && cx == dx){
      r /= 2;
    }
    if (ay == by && by == cy && cy == dy){
      g /= 2;
    }
    if (az == bz && bz == cz && cz == dz){
      b /= 2;
    }
    if ((int)z % 2){
      r /= 2;
      g /= 2;
      b /= 2;
    }
  }else if (shader == 1){
    if (ax == bx && bx == cx && cx == dx){
      r /= 2; g /= 2; b /= 2;
    }
    if (ay == by && by == cy && cy == dy){
      r /= 4; g /= 4; b /= 4;
    }
  }
  glColor3ub(r,g,b);
  glVertex3f(x+ax,y+ay,z+az);
  glColor3ub(r,g,b);
  glVertex3f(x+bx,y+by,z+bz);
  glColor3ub(r,g,b);
  glVertex3f(x+cx,y+cy,z+cz);
  glColor3ub(r,g,b);
  glVertex3f(x+dx,y+dy,z+dz);
}

void draw_vox(){
  glBegin( GL_QUADS );
  for (int i = 0; i < slice; i++){
    for (int j = 0; j < vox.h; j++){
      for (int k = 0; k < vox.w; k++){
        int a0 = rast_get(&vox,k,j,i);
        int a1 = rast_get(&vox,k-1,j,i);
        int a2 = rast_get(&vox,k+1,j,i);
        int a3 = rast_get(&vox,k,j-1,i);
        int a4 = rast_get(&vox,k,j+1,i);
        int a5 = rast_get(&vox,k,j,i-1);
        int a6 = rast_get(&vox,k,j,i+1);
        if (!a0) continue;
        if (!a1) draw_face(k,j,i,a0,0,0,1, 0,1,1, 0,1,0, 0,0,0);
        if (!a2) draw_face(k,j,i,a0,1,0,0, 1,1,0, 1,1,1, 1,0,1);
        if (!a3) draw_face(k,j,i,a0,1,0,0, 1,0,1, 0,0,1, 0,0,0);
        if (!a4) draw_face(k,j,i,a0,0,1,0, 0,1,1, 1,1,1, 1,1,0);
        if (!a6) draw_face(k,j,i,a0,1,0,1, 1,1,1, 0,1,1, 0,0,1);
        if (!a5) draw_face(k,j,i,a0,0,0,0, 0,1,0, 1,1,0, 1,0,0);
      }
    }
  }
  glEnd();
}

void draw_string(const char* text, int x, int y){
  int i = 0;
  int x0 = x;
  y -= 13;
  while (text[i]){
    glRasterPos2f(x,y);
    if (text[i] == '\n'){
      y -= 13;
      x = x0;
    }else{
      glutBitmapCharacter(GLUT_BITMAP_8_BY_13,text[i]);
      x += 8;
    }
    i++;
  }
}

void draw(){

  //select clearing (background) color
  glClearColor(0,0,0, 0.0);
  

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(60, (float)W/(float)H, 0.1, 10000);
  glEnable(GL_CULL_FACE);
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  
  float mat0[] = MAT_ROTX(rotx);
  float mat1[] = MAT_ROTY(roty);
  float mat2[] = MAT_ROTZ(rotz);
  float mat3[] = MAT_TRSL(panx,pany,panz);

  float mat5[] = MAT_MULT(mat0,mat1);
  float mat6[] = MAT_MULT(mat5,mat2);
  float mat7[] = MAT_MULT(mat6,mat3);

  glMatrixMode(GL_MODELVIEW);

  glLoadMatrixf(mat7);
  
  glEnable(GL_DEPTH_TEST);

  glCallList(dsp);

  glLoadIdentity();

  glColor3f(1.0,1.0,1.0); 
  glDisable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0,W,0,H);

  // glLoadIdentity();
  draw_string("[WASDZX] to pan\n[QERF] to rotate\n[0] to reset\n[K] to switch shader\n[M] to switch to layer view\n",2,H);

  // printf("%d\n",frame);
  
}


void draw2(){
  glClearColor(0,0,0,0.0);
  glViewport(0,0,256,256);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(60, 1, 0.1, 10000);
  glEnable(GL_CULL_FACE);
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  
  float mat0[] = MAT_ROTX(rotx);
  float mat1[] = MAT_ROTY(roty);
  float mat2[] = MAT_ROTZ(rotz);
  float mat3[] = MAT_TRSL(panx,pany,panz);

  float mat5[] = MAT_MULT(mat0,mat1);
  float mat6[] = MAT_MULT(mat5,mat2);
  float mat7[] = MAT_MULT(mat6,mat3);

  glMatrixMode(GL_MODELVIEW);

  glLoadMatrixf(mat7);
  
  glEnable(GL_DEPTH_TEST);

  glCallList(dsp);

  glBegin(GL_LINE_STRIP);
  glVertex3f(0,0,curr_slice);
  glVertex3f(vox.w,0,curr_slice);
  glVertex3f(vox.w,vox.h,curr_slice);
  glVertex3f(0,vox.h,curr_slice);
  glVertex3f(0,0,curr_slice);
  glEnd();


  glViewport(0,0,W,H);

  glLoadIdentity();

  glColor3f(1.0,1.0,1.0); 
  glDisable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0,W,0,H);
  // glLoadIdentity();

  glPushMatrix();
  
  glTranslatef(W/2,H/2,0);

  glScalef(pix_scale,pix_scale,1);

  glTranslatef(-vox.w/2,-vox.h/2,0);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  
  // glBlendFunc(GL_ZERO, GL_ONE);  

  glBegin(GL_QUADS);

  for (int i = 0; i < vox.h; i++){
    for (int j = 0; j < vox.w; j++){
      int v = rast_get(&vox,j,i,curr_slice);
      if (v){
        int rgb = colortable[(v+1) % (sizeof(colortable)/sizeof(int))];
        unsigned char r = ((rgb >> 8)  & 0xf)*16;
        unsigned char g = ((rgb >> 4)  & 0xf)*16;
        unsigned char b = ((rgb >> 0)  & 0xf)*16;
        glColor3ub(r,g,b);
      }else{
        glColor4f(0.15,0.15,0.15,0.7);
      }

      glVertex2f(j,i);
      glVertex2f(j+1,i);
      glVertex2f(j+1,i+1);
      glVertex2f(j,i+1);
    }
  }

  {
    int mx = ((mouse_x)-W/2)/pix_scale+vox.w/2;
    int my = ((H-mouse_y)-H/2)/pix_scale+vox.h/2;
    glColor4f(1.0,1.0,0.0,0.5);
    glVertex2f(mx,my);
    glVertex2f(mx+1,my);
    glVertex2f(mx+1,my+1);
    glVertex2f(mx,my+1);
    if (mouse_is_down){
      rast_set(&vox,mx,my,curr_slice,255);
    }
  }
  glEnd();
  glPopMatrix();

  glColor3f(1.0,1.0,1.0); 

  draw_string("[CV/SHIFT+CV] to iterate layers\n[L] to toggle scale\n[M] to switch to model view\n",2,H);
  

  // printf("%d\n",frame);
}

void onkeypress(unsigned char key, int x, int y){
  // printf("%d %d %d\n",(int)key,x,y);

  if (key == 27){
    exit(0);
  }else if (key == 's'){
    panz -= 100;
  }else if (key == 'w'){
    panz += 100;
  }else if (key == 'a'){
    panx += 100;
  }else if (key == 'd'){
    panx -= 100;
  }else if (key == 'z'){
    pany -= 100;
  }else if (key == 'x'){
    pany += 100;
  }else if (key == 'q'){
    roty -= 0.1;
  }else if (key == 'e'){
    roty += 0.1;
  }else if (key == 'r'){
    rotx += 0.1;
  }else if (key == 'f'){
    rotx -= 0.1;
  }else if (key == '0'){
    rotx = -M_PI/2;
    roty = 0;
    rotz = 0;
    panx = -vox.w/2;
    panz = -vox.h/2-fmax(vox.w,vox.d);
    pany = -vox.d/2;
  }
  
  // if (key == 'c'){
  //   slice = MAX(0,slice-1);
  //   glNewList(dsp, GL_COMPILE);
  //   draw_vox();
  //   glEndList();
  // }else if (key == 'v'){
  //   slice = MIN(vox.d,slice+1);
  //   glNewList(dsp, GL_COMPILE);
  //   draw_vox();
  //   glEndList();
  // }

  if (key == 'c'){
    curr_slice = MAX(0,curr_slice-1);
    printf("z=%d\n",curr_slice);
  }else if (key == 'v'){
    curr_slice = MIN(vox.d,curr_slice+1);
    printf("z=%d\n",curr_slice);
  }else if (key == 'C'){
    curr_slice = MAX(0,curr_slice-20);
    printf("z=%d\n",curr_slice);
  }else if (key == 'V'){
    curr_slice = MIN(vox.d,curr_slice+20);
    printf("z=%d\n",curr_slice);
  }else if (key == 'm'){
    mode = !mode;
  }else if (key == 'l'){
    pix_scale = 3-pix_scale;
  }

  if (key == 'k'){
    shader = !shader;
    printf("[display list] regenerating...\n");
    glNewList(dsp, GL_COMPILE);
    draw_vox();
    glEndList();
    printf("[display list] regenerated\n");
  }
}

void onmouseclick(int btn, int state, int x, int y){
  printf("%d %d %d %d\n",btn, state, x,y);
  mouse_is_down = !state;
  mouse_x = x;
  mouse_y = y;
}

void onmousemove(int x, int y){
  mouse_is_down = 0;
  mouse_x = x;
  mouse_y = y;
}

void onmousedrag(int x, int y){
  mouse_is_down = 1;
  mouse_x = x;
  mouse_y = y;
}

void display(void){
  if (mode == 0){
    draw();
  }else{
    draw2();
  }
  glFlush();
  frame++;
}
void animation_frame(){
  glutPostRedisplay();
  glutTimerFunc( 10, animation_frame, 1);
}
void onreshape(){
  W = glutGet(GLUT_WINDOW_WIDTH);
  H = glutGet(GLUT_WINDOW_HEIGHT);

}


int main(int argc, char** argv){
  if (argc < 2){
    printf("usage: voxel_viewer input.bin\n");
    printf("( voxel binary format: (LE)                  \n");
    printf("w (u32), h (u32), d (u32), data (u8,u8,u8,...\n");
    printf("voxel(x,y,z) = data[z*w*h + y*w + x] )       \n");
    exit(0);
  }

  vox = read_voxels_from_file(argv[1]);

  rotx = -M_PI/2;
  panx = -vox.w/2;
  panz = -vox.h/2-fmax(vox.w,vox.d);
  pany = -vox.d/2;
  slice = vox.d;

  if (argc > 2){
    slice = atoi(argv[2]);
  }

  glutInit(&argc, argv);  
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(W,H);
  glutInitWindowPosition(0,0);
  glutCreateWindow(" ");
  
  printf("[display list] generating...\n");

  dsp = glGenLists (1);

  glNewList(dsp, GL_COMPILE);
  draw_vox();
  glEndList();

  printf("[display list] generated.\n");

  glutDisplayFunc(display);
  glutReshapeFunc(onreshape);
  glutKeyboardFunc(onkeypress);
  glutMouseFunc(onmouseclick);
  glutPassiveMotionFunc(onmousemove);
  glutMotionFunc(onmousedrag);
  glutTimerFunc( 1, animation_frame, 1);
  
  glutMainLoop();
  return 0;
}

