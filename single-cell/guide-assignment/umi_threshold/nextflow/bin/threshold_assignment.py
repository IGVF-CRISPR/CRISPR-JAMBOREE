import argparse
import unittest

import mudata as md
from scipy.sparse import dok_matrix


def threshold_posteriors(guides, threshold):
    guide_count_array = guides.X.todok()
    threshold = 5 if threshold is None else threshold
    array = dok_matrix(guides.X.shape)
    for (x, y), guide_count in guide_count_array.items():
        if guide_count >= threshold:
            array[x, y] = 1
    return array.tocsr()


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", required=True, type=str, help="Input MuData file"
    )
    parser.add_argument(
        "-o", "--output", required=True, type=str, help="Output MuData file"
    )
    parser.add_argument("-t", "--threshold", default=None, type=float)

    return parser.parse_args()


def run(gas, threshold):
    guides = gas["guide"]
    guides.layers["guide_assignment"] = threshold_posteriors(guides, threshold)


if __name__ == "__main__":
    args = get_args()
    mu_input = md.read(args.input)
    run(mu_input, args.threshold)
    md.write(args.output, mu_input)


####################################
#              Tests               #
####################################


class TestThreshold(unittest.TestCase):
    """Test guide assignment using basic thresholding"""

    def test_layer_values(self):
        """Ensure all the values in the guide_assignment layer are 1.0. It's a sparse matrix so any locations without values
        are implicitly 0.0"""

        gas = md.read("test_data/gasperini_guide_assignment_input_minimal.h5mu")
        run(gas, threshold=5)
        guides = gas["guide"]
        self.assertIn("guide_assignment", guides.layers)
        assignments = guides.layers[
            "guide_assignment"
        ].tocoo()  # Must be in COO format to iterate overvalues
        for value in assignments.data:
            self.assertEqual(value, 1.0)
