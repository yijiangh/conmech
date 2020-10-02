import pytest
import os, glob
from pyconmech import StiffnessChecker

@pytest.mark.skip(reason='not fully developed')
@pytest.mark.subd
def test_subdivided_compliance():
    here = os.path.dirname(os.path.abspath(__file__))
    filenames_num = glob.glob(os.path.join(here, '..', 'test_data', 'subdivided', 'topopt100_subd_num'+'*'+'.json'))
    filenames_len = glob.glob(os.path.join(here, '..', 'test_data', 'subdivided', 'topopt100_subd_num'+'*'+'.json'))
    filenames = filenames_num + filenames_len

    for file in filenames:
        print('\n'+ '*' * 10)
        print('test file: {}'.format(file))
        sc = StiffnessChecker.from_json(json_file_path=file)
        print('num of pts: {}, num of elements: {}'.format(len(sc.node_points), len(sc.elements)))
        sc.set_loads(gravity_direction=[0,0,-1.0])
        success = sc.solve()
        print('success: {}'.format(success))
        print('has stored result: {}'.format(sc.has_stored_result()))
        print('compliance: {}'.format(sc.get_compliance()))


if __name__ == '__main__':
    test_subdivided_compliance()