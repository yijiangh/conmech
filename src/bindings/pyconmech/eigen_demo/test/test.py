import eigen_demo_py as m

p = m.EigenSolveDemo()

print(p)

json_path =\
    "/home/yijiangh/Documents/assembly-instances/assembly_models/spatial_extrusion/voronoi/voronoi_S1.0_09-05-2018.json"

p.testEigen()
p.testJsonParse(json_path)