import cgamma
import numpy


class Parameters:
    def __init__(self,
                 diff:      float = 0.03,
                 dta:       float = 2.0,
                 threshold: float = 0.10,
                 norm:      str   = "GLOBAL",
                 relative:  bool  = False):
        self.diff = diff
        self.dta = dta
        self.threshold = threshold
        self.norm = norm
        self.relative = relative


class Options:
    def __init__(self, pass_only: bool = False, pattern_shrinks: int = 6):
        self.pass_only = pass_only
        self.pattern_shrinks = pattern_shrinks


class Distribution:
    def __init__(self):
        self.matrix = numpy.identity(3, dtype=numpy.double)
        self.origin = numpy.zeros((3, 1), dtype=numpy.double)
        self.spacing = numpy.ones((3, 1), dtype=numpy.double)
        self.data: numpy.ndarray = None


class Results:
    def __init__(self):
        self.total = 0
        self.passed = 0
        self.min = 0.0
        self.max = 0.0
        self.mean = 0.0
        self.msqr = 0.0
        self.dist: numpy.ndarray = None


def compute(params:  Parameters,
            options: Options,
            ref:     Distribution,
            meas:    Distribution):
    res = Results()
    cgamma.compute(params, options, ref, meas, res)
    return res
