"""Load the full 21-study NIDM pain studyset with images."""
import os, warnings
import nimare
from nimare.transforms import ImageTransformer
from nimare.tests.utils import get_test_data_path
warnings.simplefilter("ignore")

def load_pain():
    dset = nimare.dataset.Dataset(
        os.path.join(get_test_data_path(), "nidm_pain_dset.json"))
    dset.update_path("/root/.nimare/nidm_21pain")
    ss = nimare.studyset.normalize_collection(dset)
    return ImageTransformer(target="z").transform(ss)
