all: lucy bunny nefertiti dragon knot23 knot35 boy

folders:
	mkdir -p models
	mkdir -p output

clean:
	rm -rf models
	rm -rf output

lucy: folders
	[ -f models/lucy.obj ] || curl -L "https://cdn.glitch.com/4c1e3290-9021-4225-bce9-a7df311fcec1%2Flucy_100K.obj?v=1601077747141" > models/lucy.obj
	../mesh_to_voxels models/lucy.obj -o output/lucy.bin --inv_flood 1 --swap_yz
	../accentuate_voxels output/lucy.bin -o output/lucy-acc.bin --zero_sum
	../voxels_to_gcode output/lucy-acc.bin -o output/lucy-acc.gcode --config ../profiles/ender3.ini --sink 8 --step 5

bunny: folders
	[ -f models/bunny.obj ] || curl -L "https://cdn.glitch.global/4c1e3290-9021-4225-bce9-a7df311fcec1/bunny_70K.obj?v=1663891041460" > models/bunny.obj
	../mesh_to_voxels models/bunny.obj -o output/bunny.bin --inv_flood 1 --step 0.1 --swap_yz --fit_z 100
	../accentuate_voxels output/bunny.bin -o output/bunny-acc.bin --zero_sum --smooth 6 --rate 9
	../voxels_to_gcode output/bunny-acc.bin -o output/bunny-acc.gcode --config ../profiles/ender3.ini --step 5 --sink 15

nefertiti: folders
	[ -f models/nefertiti.obj ] || curl -L "https://cdn.glitch.global/4c1e3290-9021-4225-bce9-a7df311fcec1/nefertiti_invert_100K.obj?v=1663892024244" > models/nefertiti.obj
	../mesh_to_voxels models/nefertiti.obj -o output/nefertiti.bin --inv_flood 2 --fit_z 100
	../accentuate_voxels output/nefertiti.bin -o output/nefertiti-acc.bin --rate 5 --smooth 7 --zero_sum
	../voxels_to_gcode output/nefertiti-acc.bin -o output/nefertiti-acc.gcode --config ../profiles/ender3.ini --step 4 --sink 8

dragon: folders
	[ -f models/dragon.obj ] || curl -L "https://cdn.glitch.com/4c1e3290-9021-4225-bce9-a7df311fcec1%2Fdragon_100K.obj?v=1601077751925" > models/dragon.obj
	../mesh_to_voxels models/dragon.obj -o output/dragon.bin --inv_flood 3 --fit_z 68 --swap_yz
	../accentuate_voxels output/dragon.bin -o output/dragon-acc.bin --rate 7 --zero_sum
	../generate_support output/dragon-acc.bin -o output/dragon-acc-sup.bin --sink 8 --hang 5 --ero_xy 5
	../voxels_to_gcode output/dragon-acc-sup.bin -o output/dragon-acc-sup.gcode --config ../profiles/ender3.ini --step 3 --sink 18

knot23: folders
	gcc -O3 knot.c -o output/knot
	./output/knot 2 3
	../voxels_to_gcode output/knot23.bin -o output/knot23.gcode --config ../profiles/ender3.ini --step 4 --sink 4

knot35: folders
	gcc -O3 knot.c -o output/knot
	./output/knot 3 5
	../voxels_to_gcode output/knot35.bin -o output/knot35.gcode --config ../profiles/ender3.ini --step 4 --sink 4

boy: folders
	gcc -O3 boy.c -o output/boy
	./output/boy
	../voxels_to_gcode output/boy-solid.bin -o output/boy-solid.gcode --config ../profiles/ender3.ini --step 4 --sink 2