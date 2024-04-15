"""folders.py
Variables to store the folder names
"""

RESULTS = Path("results")

# RAW = RESULTS / "raw"

READS = RESULTS / "reads"
REFERENCE = RESULTS / "reference"

PRE = RESULTS / "preprocess"
PRE_INDEX = PRE / "index"
PRE_MAP = PRE / "map"
PRE_SPLIT = PRE / "split"
PRE_FILT = PRE / "filtered"
PRE_MPILEUP = PRE / "mpileup"

# mpileup - SNP calling - popoolation
POP1 = RESULTS / "popoolation1"
POP1_RAW = POP1 / "raw"
POP1_FILT = POP1 / "filtered"
POP1_SUB = POP1 / "subsampled"
POP1_TABLES = POP1 / "tables"
POP1_PLOTS = POP1 / "plots"

# sync - popoolation2
POP2 = RESULTS / "popoolation2"
POP2_MPILEUP = POP2 / "mpileup"
POP2_FILT = POP2 / "filt"
POP2_SYNC = POP2 / "sync"  # rename?
POP2_SUB = POP2 / "subsampled"
POP2_TABLES = POP2 / "tables"
POP2_PLOTS = POP2 / "plots"

# Tables & plots
# POPOOLATION = RESULTS / "popoolation"
# POPOOLATION_TABLES = POPOOLATION / "tables"
# POPOOLATION_PLOTS = POPOOLATION / "plots"

HP = RESULTS / "hp"
HP_TABLES = HP / "tables"
HP_PLOTS = HP / "plots"
# FST = RESULTS / "fst"
# FST_TABLES = FST / "tables"
# FST_PLOTS = FST / "plots"
