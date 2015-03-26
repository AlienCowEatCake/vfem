Merge "cube.geo";
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh 3;
Save Sprintf("cube_x%g.msh", 1);
Mesh.SecondOrderLinear = 1;
For i In {2:4:2} 
  RefineMesh; 
  Save Sprintf("cube_x%g.msh", i); 
EndFor
