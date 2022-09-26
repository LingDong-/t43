# -Wall -std=c99 -Wno-unused-variable
all:
	gcc -O3 mesh_to_voxels.c -o mesh_to_voxels
	gcc -O3 march_voxels.c -o march_voxels
	gcc -O3 blur_voxels.c -o blur_voxels
	gcc -O3 voxels_to_gcode.c -o voxels_to_gcode
	gcc -O3 generate_support.c -o generate_support
	gcc -O3 voxel_convert.c -o voxel_convert
	gcc -O3 accentuate_voxels.c -o accentuate_voxels

viewer:
	gcc -O3 voxel_viewer.c -framework OpenGL -framework GLUT -o voxel_viewer
  #gcc -O3 voxel_viewer.c -lgl -lglut -o voxel_viewer