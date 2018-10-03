import conmech_py as cm

json_path =\
    "/home/yijiangh/Documents/assembly-instances/assembly_models/spatial_extrusion/voronoi/voronoi_S1.0_09-05-2018.json"


sc = cm.stiffness_checker(json_path, True)
ids = [1,2,3]
sc.check_deformation(ids)
