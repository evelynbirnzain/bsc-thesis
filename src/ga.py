from tqdm import tqdm

import numpy as np
import pandas as pd
import os
import io_utils
from scoring import Scoring
from plotting import Plotting


class GA:


    def __init__(self, genes: [str], scoring: Scoring, plotting: Plotting, init_pop_weights=None):
        self.genes = genes
        self.init_range_high = len(self.genes) - 1
        self.gene_space = range(0, len(self.genes))
        self.random_mutation_max_val = len(self.genes) - 1
        self.plotting = plotting
        self.scoring = scoring
        self.init_pop_weights = init_pop_weights








