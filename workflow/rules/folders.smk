"""folders.py
Variables to store the folder names
"""

RESULTS = Path("results")

# Raw data
READS = RESULTS / "reads"
REFERENCE = RESULTS / "reference"

# preprocess
PRE = RESULTS / "preprocess"
PRE_INDEX = PRE / "index"
PRE_MAP = PRE / "map"
PRE_SPLIT = PRE / "split"
PRE_FILT = PRE / "filtered"
PRE_MPILEUP = PRE / "mpileup"

# popoolation - sample-wise
POP1 = RESULTS / "popoolation1"
POP1_RAW = POP1 / "raw"
POP1_FILT = POP1 / "filtered"
POP1_SUB = POP1 / "subsampled"
POP1_VS = POP1 / "variance_sliding"
POP1_PLOTS = POP1 / "plots"
POP1_HP = POP1 / "hp"

# popoolation2 - pair-wise
POP2 = RESULTS / "popoolation2"
POP2_MPILEUP = POP2 / "mpileup"
POP2_FILT = POP2 / "filtered"
POP2_SYNC = POP2 / "sync"
POP2_SUB = POP2 / "subsampled"
POP2_FST = POP2 / "fst_sliding"
POP2_PLOTS = POP2 / "plots"
