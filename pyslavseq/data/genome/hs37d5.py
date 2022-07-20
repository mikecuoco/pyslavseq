#!/usr/bin/env python
__author__ = "Michael Cuoco"

from ...genome import Genome
from .. import genome as _myself

hs37d5 = Genome(_myself.__path__[0] + '/hs37d5.genome')
