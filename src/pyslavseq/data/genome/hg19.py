#!/usr/bin/env python
__author__ = "Apuã Paquola"

from ...genome import Genome
from .. import genome as _myself

hg19 = Genome(_myself.__path__[0] + '/hg19.genome')
