import numpy as np
from termcolor import cprint
from pyconmech.frame_analysis.stiffness_checker import StiffnessChecker
from pyconmech.frame_analysis.io_base import Model, LoadCase

def solve_linear_elastic_from_data(model_data, loadcase_data, verbose=True):
    model = Model.from_data(model_data)
    loadcase = LoadCase.from_data(loadcase_data)

    sc = StiffnessChecker(model, checker_engine="numpy", verbose=verbose)

    sc.set_loads(loadcase)
    
    # if the structure's nodal deformation exceeds 1e-3 meter,
    # we want the checker to return `sol_success = False`
    sc.set_nodal_displacement_tol(trans_tol=np.inf, rot_tol=np.inf)
    
    success = sc.solve()

    # Get all the analysis information:
    # nodal deformation, fixity reactions, element reactions
    success, nD, fR, eR = sc.get_solved_results()
    cprint('Solve success: {}'.format(success), 'green' if success else 'red')

    # nD is the nodal displacements
    # a dictionary indexed by nodal indices, values are the [ux, uy, uz, theta_x, theta_y, theta_z] deformation vector
    trans_tol, rot_tol = sc.get_nodal_deformation_tol()
    max_trans, max_rot, max_trans_vid, max_rot_vid = sc.get_max_nodal_deformation()
    # translation = np.max(np.linalg.norm([d[:3] for d in nD.values()], ord=2, axis=1))
    # rotation = np.max(np.linalg.norm([d[3:] for d in nD.values()], ord=2, axis=1))
    cprint('Max translation deformation: {0:.5f} m / {1:.5} = {2:.5}, at node #{3}'.format(
        max_trans, trans_tol, max_trans / trans_tol, max_trans_vid), 'cyan')
    cprint('Max rotation deformation: {0:.5f} rad / {1:.5} = {2:.5}, at node #{3}'.format(
        max_rot, rot_tol, max_rot / rot_tol, max_rot_vid), 'cyan')
    
    # Compliance is the elastic energy: https://en.wikipedia.org/wiki/Elastic_energy
    # The inverse of stiffness is flexibility or compliance
    # the higher this value is, the more flexible a structure is
    # we want this value to be low
    compliance = sc.get_compliance()
    cprint('Elastic energy: {}'.format(compliance), 'cyan')
    # volume can be computed by simply summing over `cross sectional area * bar length`

    rdata = {
        'solve_success' : bool(success),
        'elastic_energy' : float(compliance),
        'nodal_displacement' : { nid : list(nd) for nid, nd in nD.items()},
        'support_reaction' :   { nid : list(fr) for nid, fr in fR.items()},
        'element_reaction' :   { eid : [list(er[0]), list(er[1])] for eid, er in eR.items()},
        'max_trans' : max_trans,
        'max_trans_nid' : int(max_trans_vid),
    }

    # if write:
    #     result_path = os.path.abspath(os.path.join(problem_path, problem + '_result.json'))
    #     with open(result_path, 'w') as f:
    #         json.dump(rdata, f, indent=None)
    #     cprint('Analysis result saved to {}'.format(result_path), 'green')
    return rdata
    