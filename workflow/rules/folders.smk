"""folders.py
Variables to store the folder names
"""

RESULTS = Path("results")

# RAW = RESULTS / "raw"

READS = RESULTS / "reads"
REFERENCE = RESULTS / "reference"

PRE = RESULTS / "preprocess"
QC = PRE / "qc"
MAP = PRE / "map"
MAP_INDEX = MAP / "index"
MAP_RAW = MAP / "raw"
MAP_RAW = MAP / "raw"
MAP_SPLIT = MAP / "split"
MAP_FILT = MAP / "filtered"
PRE_MPILEUP = PRE / "mpileup"

# mpileup - SNP calling - popoolation
POP1 = RESULTS / "popoolation1"
MPILEUP = RESULTS / "mpileup"
MPILEUP_RAW = MPILEUP / "raw"
MPILEUP_FILT = MPILEUP / "filt"
MPILEUP_SUB = MPILEUP / "sub"

# sync - popoolation2
POP2 = RESULTS / "popoolation2"
SYNC = RESULTS / "sync"
SYNC_MPILEUP = SYNC / "mpileup"
SYNC_FILT = SYNC / "filt"
SYNC_MPILEUP2SYNC = SYNC / "mpileup2sync"  # rename?
SYNC_SUBSAMPLED = SYNC / "subsampled"



# Tables & plots
POPOOLATION = RESULTS / "popoolation"
POPOOLATION_TABLES = POPOOLATION / "tables"
POPOOLATION_PLOTS = POPOOLATION / "plots"

HP = RESULTS / "hp"
HP_TABLES = HP / "tables"
HP_PLOTS = HP / "plots"

FST = RESULTS / "fst"
FST_TABLES = FST / "tables"
FST_PLOTS = FST / "plots"
