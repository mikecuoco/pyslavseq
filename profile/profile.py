import poetry
import cProfile
import pstats
from pstats import SortKey

pr = cProfile.Profile()
pr.enable()

# pyslavseq_extract_features {arguments}

pr.disable()
stats = pstats.Stats(pr).sort_stats('cumtime')
stats.dump_stats("profile/results")