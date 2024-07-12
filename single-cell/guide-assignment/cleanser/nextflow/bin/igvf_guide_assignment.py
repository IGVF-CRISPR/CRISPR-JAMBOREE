import argparse
import asyncio
import unittest

import mudata as md
import numpy as np
from cleanser import CS_MODEL_FILE, DC_MODEL_FILE
from cleanser import run as run_cleanser
from scipy.sparse import dok_matrix


def posteriors_layer(stan_results, array, threshold=None):
    if threshold is None:
        for guide_id, (samples, cell_info) in stan_results.items():
            pzi = np.transpose(samples.stan_variable("PZi"))
            for i, (cell_id, _) in enumerate(cell_info):
                array[cell_id, guide_id] = np.median(pzi[i])
    else:
        for guide_id, (samples, cell_info) in stan_results.items():
            pzi = np.transpose(samples.stan_variable("PZi"))
            for i, (cell_id, _) in enumerate(cell_info):
                if np.median(pzi[i]) >= threshold:
                    array[cell_id, guide_id] = 1

    return array.tocsr()


def cleanser_posteriors(guides, threshold):
    guide_count_array = guides.X.todok()
    counts = [
        (key[1], key[0], int(guide_count))
        for key, guide_count in guide_count_array.items()
    ]
    analysis = guides.uns.get("capture_method")
    if analysis is None or analysis[0] == "CROP-seq":
        model = CS_MODEL_FILE
    elif analysis == "direct capture":
        model = DC_MODEL_FILE
    else:
        raise ValueError("Invalid capture method type")

    results = asyncio.run(run_cleanser(counts, model))
    return posteriors_layer(results, dok_matrix(guides.X.shape), threshold)


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

    model_group = parser.add_mutually_exclusive_group(required=True)
    model_group.add_argument(
        "--cleanser",
        action="store_true",
        help="Use CLEANSER to determine assignments",
    )
    model_group.add_argument(
        "--umi-threshold",
        action="store_true",
        help="Use UMI threshold to determine assignments (default 5)",
    )
    return parser.parse_args()


def run(gas, cleanser, umi_threshold, threshold):
    guides = gas["guide"]

    if cleanser:
        posteriors = cleanser_posteriors(guides, threshold)
    elif umi_threshold:
        posteriors = threshold_posteriors(guides, threshold)

    guides.layers["guide_assignment"] = posteriors


if __name__ == "__main__":
    args = get_args()
    mu_input = md.read(args.input)
    run(mu_input, args.cleanser, args.umi_threshold, args.threshold)
    md.write(args.output, mu_input)


####################################
#              Tests               #
####################################


class TestCleanser(unittest.TestCase):
    """Test guide assignment using basic thresholding"""

    def test_cropseq_threshold_values(self):
        """Ensure all the values in the guide_assignment layer are 1.0. It's a sparse matrix so any locations without values
        are implicitly 0.0"""

        gas = md.read("test_data/gasperini_guide_assignment_input_minimal.h5mu")
        guides = gas["guide"]

        self.assertEqual("CROP-seq", guides.uns.get("capture_method"))

        run(gas, cleanser=True, umi_threshold=False, threshold=0.8)
        self.assertIn("guide_assignment", guides.layers)
        assignments = guides.layers[
            "guide_assignment"
        ].tocoo()  # Must be in COO format to iterate overvalues
        for value in assignments.data:
            self.assertEqual(value, 1.0)

    def test_cropseq_probability_values(self):
        """Ensure all the values in the guide_assignment layer are > 0.0. It's a sparse matrix so any locations without values
        are implicitly 0.0"""

        gas = md.read("test_data/gasperini_guide_assignment_input_minimal.h5mu")
        guides = gas["guide"]

        self.assertEqual("CROP-seq", guides.uns.get("capture_method"))

        run(gas, cleanser=True, umi_threshold=False, threshold=None)
        self.assertIn("guide_assignment", guides.layers)
        assignments = guides.layers[
            "guide_assignment"
        ].tocoo()  # Must be in COO format to iterate overvalues

        # All values are > 0
        for value in assignments.data:
            self.assertGreater(value, 0.0)

        # Not all values are == 1.0
        self.assertFalse(all([d == 1.0 for d in assignments.data]))

    def test_dc_threshold_values(self):
        """Ensure all the values in the guide_assignment layer are 1.0. It's a sparse matrix so any locations without values
        are implicitly 0.0"""

        gas = md.read("test_data/papalexi_guide_assignment_input.h5mu")
        guides = gas["guide"]

        self.assertEqual("direct capture", guides.uns.get("capture_method"))

        run(gas, cleanser=True, umi_threshold=False, threshold=0.8)
        self.assertIn("guide_assignment", guides.layers)
        assignments = guides.layers[
            "guide_assignment"
        ].tocoo()  # Must be in COO format to iterate overvalues
        for value in assignments.data:
            self.assertEqual(value, 1.0)

    def test_dc_probability_values(self):
        """Ensure all the values in the guide_assignment layer are > 0.0. It's a sparse matrix so any locations without values
        are implicitly 0.0"""

        gas = md.read("test_data/papalexi_guide_assignment_input.h5mu")
        guides = gas["guide"]

        self.assertEqual("direct capture", guides.uns.get("capture_method"))

        run(gas, cleanser=True, umi_threshold=False, threshold=None)
        self.assertIn("guide_assignment", guides.layers)
        assignments = guides.layers[
            "guide_assignment"
        ].tocoo()  # Must be in COO format to iterate overvalues

        # All values are > 0
        for value in assignments.data:
            self.assertGreater(value, 0.0)

        # Not all values are == 1.0
        self.assertFalse(all([d == 1.0 for d in assignments.data]))


class TestThreshold(unittest.TestCase):
    """Test guide assignment using basic thresholding"""

    def test_layer_values(self):
        """Ensure all the values in the guide_assignment layer are 1.0. It's a sparse matrix so any locations without values
        are implicitly 0.0"""

        gas = md.read("test_data/gasperini_guide_assignment_input_minimal.h5mu")
        run(gas, cleanser=False, umi_threshold=True, threshold=5)
        guides = gas["guide"]
        self.assertIn("guide_assignment", guides.layers)
        assignments = guides.layers[
            "guide_assignment"
        ].tocoo()  # Must be in COO format to iterate overvalues
        for value in assignments.data:
            self.assertEqual(value, 1.0)
