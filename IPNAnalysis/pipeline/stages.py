# ==========================================================
# pipeline/stages.py
# ProcessingStage enum, NpEncoder, and checkpoint schema constants.
# ==========================================================

import json
from enum import Enum

import numpy as np


class NpEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy scalar and array types."""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


class ProcessingStage(Enum):
    NOT_STARTED = 0
    PREPROCESSING = 1
    PREPROCESSING_COMPLETE = 2
    SORTING = 3
    SORTING_COMPLETE = 4
    MERGE = 5
    MERGE_COMPLETE = 6
    ANALYZER = 7
    ANALYZER_COMPLETE = 8
    REPORTS = 9
    REPORTS_COMPLETE = 10


#: Increment this when the checkpoint schema changes so old checkpoints can be
#: migrated automatically on load.
CHECKPOINT_SCHEMA_VERSION = 2
