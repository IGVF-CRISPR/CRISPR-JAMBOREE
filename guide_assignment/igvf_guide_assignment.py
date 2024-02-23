import argparse
import asyncio

import mudata as md
import numpy as np
from scipy.sparse import dok_matrix

from cleanser import CS_MODEL_FILE, DC_MODEL_FILE, run


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


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="Input MuData file")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output MuData file")
    parser.add_argument("-t", "--threshold", default=None, type=float)
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    gas = md.read(args.input)
    guides = gas["guide"]
    count_array = guides.X.todok()
    counts = [(key[1], key[0], int(value)) for key, value in count_array.items()]

    analysis = guides.uns.get("capture_method")
    if analysis is None or analysis[0] == "CROP-seq":
        model = CS_MODEL_FILE
    elif analysis == "direct capture":
        model = DC_MODEL_FILE
    else:
        raise ValueError("Invalid capture method type")

    results = asyncio.run(run(counts, model))
    posteriors = posteriors_layer(results, dok_matrix(guides.X.shape), args.threshold)
    guides.layers["guide_assignment"] = posteriors
    md.write(args.output, gas)
