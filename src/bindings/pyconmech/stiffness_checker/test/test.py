import conmech_py as cm

json_path =\
    "/home/yijiangh/Documents/assembly-instances/assembly_models/spatial_extrusion/voronoi/voronoi_S1_09-26-2018.json"


sc = cm.stiffness_checker(json_path, False)

ids = [1,2,3]
sc.check_deformation(ids)