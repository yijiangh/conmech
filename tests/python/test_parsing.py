import os
import sys
if sys.version_info[0] < 3:
    from backports import tempfile
else:
    import tempfile
from numpy.testing import assert_equal, assert_almost_equal
import pytest
import random

from pyconmech import StiffnessChecker
from pyconmech.frame_analysis import read_frame_json, write_frame_json, read_load_case_json, check_material_dict
from copy import deepcopy

@pytest.mark.parse
def test_frame_file_io(test_data_dir):
    file_name = 'tower_3D_broken_lines.json'
    file_path = os.path.join(test_data_dir, file_name)
    node_points, element_vids, fix_specs, model_type, material_dicts, model_name, unit = \
        read_frame_json(file_path, verbose=True)

    # temp_fp = os.path.join(here, 'tmp.json')
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_fp = os.path.join(temp_dir, file_name)
        write_frame_json(temp_fp, node_points, element_vids, fix_specs, material_dicts, \
            model_type=model_type, model_name=model_name, unit=unit)
        back_node_points, back_element_vids, back_fix_specs, back_model_type, back_material_dicts, back_model_name, back_unit = \
            read_frame_json(temp_fp, verbose=True)

    for n1, n2 in zip(node_points, back_node_points):
        assert_equal(n1, n2)
    for e1, e2 in zip(element_vids, back_element_vids):
        assert_equal(e1, e2)
    for vid, spec in fix_specs.items():
        assert vid in back_fix_specs
        assert_equal(spec, back_fix_specs[vid])
    assert model_type == back_model_type
    for mat1, mat2 in zip(material_dicts, back_material_dicts):
        assert mat1 == mat2
    assert model_name == back_model_name
    assert unit == back_unit

@pytest.mark.parse_mat
def test_parse_material_properties_from_frame_json(test_data_dir):
    file_name = 'bad_material_properties_model.json'
    json_path = os.path.join(test_data_dir, file_name)

    node_points, element_vids, fix_specs, model_type, material_dicts, model_name, unit = \
        read_frame_json(json_path)

    # initial file no problem
    mat_entry = ''
    err_msg = 'Missing attributes in json: uniform_cross_section or uniform_material_properties!'
    with pytest.raises(RuntimeError) as excinfo:
        StiffnessChecker.from_json(json_file_path=json_path)
    assert str(excinfo.value) == err_msg

    def check_mat(i, mat_entry, err_msg, expect_failure=True):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_fp = os.path.join(temp_dir, file_name)
            mds = deepcopy(material_dicts)
            if i is not None:
                uniform = False
                err_msg = 'Element {}: '.format(i) + err_msg
            else:
                i = 0
                uniform = True

            if isinstance(mat_entry, list):
                for m in mat_entry:
                    del mds[i][m]
            else:
                del mds[i][mat_entry]

            write_frame_json(temp_fp, node_points, element_vids, fix_specs, mds, \
              unif_cross_sec=uniform, unif_material=uniform, \
              model_type=model_type, model_name=model_name, check_material=False, unit=unit)
            if expect_failure:
                with pytest.raises(RuntimeError) as excinfo:
                    StiffnessChecker.from_json(json_file_path=temp_fp)
                # assert str(excinfo.value) == err_msg
            else:
                StiffnessChecker.from_json(json_file_path=temp_fp)

    for test_element in [True, False]:
        i = random.choice(list(range(len(element_vids)))) if test_element else None
        mat_entry = 'youngs_modulus'
        err_msg = 'Young\'s modulus property value not specified!'
        check_mat(i, mat_entry, err_msg)

        for j_name in ['density', 'density_unit']:
            i = random.choice(list(range(len(element_vids)))) if test_element else None
            mat_entry = j_name
            err_msg = 'Density not specified!'
            check_mat(i, mat_entry, err_msg)

        i = random.choice(list(range(len(element_vids)))) if test_element else None
        mat_entry = 'cross_sec_area'
        err_msg = 'Cross section area property value not specified!'
        check_mat(i, mat_entry, err_msg)

        for j_name in ['Jx', 'Iy', 'Iz']:
            i = random.choice(list(range(len(element_vids)))) if test_element else None
            mat_entry = j_name
            err_msg = 'Jx, Iy or Iz not specified!'
            check_mat(i, mat_entry, err_msg)

        # specify one of them is sufficient
        for j_name in ['poisson_ratio', 'shear_modulus']:
            i = random.choice(list(range(len(element_vids)))) if test_element else None
            mat_entry = j_name
            err_msg = ''
            check_mat(i, mat_entry, err_msg, expect_failure=False)

        i = random.choice(list(range(len(element_vids)))) if test_element else None
        mat_entry = ['shear_modulus', 'poisson_ratio']
        err_msg = 'Both shear modulus and poisson_ratio not specified!'
        check_mat(i, mat_entry, err_msg, expect_failure=True)

@pytest.mark.parse_element
def test_parse_element_from_json(test_data_dir):
    file_name = 'bad_node_element_model.json'
    json_path = os.path.join(test_data_dir, file_name)

    node_points, element_vids, fix_specs, model_type, material_dicts, model_name, unit = \
        read_frame_json(json_path)

    with pytest.raises(RuntimeError) as excinfo:
        StiffnessChecker.from_frame_data(node_points, element_vids, fix_specs, material_dicts, \
            unit='meter', model_type='frame', model_name=None, verbose=False)
    assert str(excinfo.value) == 'there needs to be at least one support (fixed) vertex in the model!'