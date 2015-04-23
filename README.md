# Mesh.Hex1
Hex1 P1 Hexahedral structured mesh rectangle (MATLAB)

    H = Hex1(x0,x1,nxe,y0,y1,nye,z0,z1,nze)

## Methods
**Hex1** includes the following methods:
- Hex1() - The constructor, creates the Hex1Mesh object.
- Neighbors() - Returns a m-by-4 matrix containing element indices to the neighbors.
- vizMesh()
- RefineLocal()
- RefineLocal()
- boundaryInds()
- baseHexP1()

## Installation
Put the contents of this repo in the *+Mesh\@Hex1* folder somewhere in your MATLAB path.

## Demo
```matlab
%Create mesh
ne = 2;
x0 = -1;x1 = 1;y0 = -1;y1 = 1;z0 = -1;z1 = 1;
nxe = ne;nye = ne;nze = ne;
H = Hex1Mesh(x0,x1,nxe,y0,y1,nye,z0,z1,nze);
H.vizMesh('ElementNumbers','NodeNumbers');
```

## Mesh numbering
The mesh is numbered according to the figure below. The red and cyan colors are used to explain the refinement algorithm.
![](http://i.imgur.com/M8ZaWMv.png)

