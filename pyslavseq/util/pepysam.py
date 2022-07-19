#!/usr/bin/env python
__author__ = "Apuã Paquola"


def pesamfile(samfile, filter_func=lambda x: True):
    """ Iterates through a paired-end sam file yielding (read1, read2) tuples """
    with samfile:
        try:
            i = 0
            previous_read = None
            while True:
                read = next(samfile)
                if not filter_func(read):
                    continue

                if i == 1:
                    assert previous_read.qname == read.qname
                    assert previous_read.is_read1
                    assert read.is_read2
                    yield (previous_read, read)

                previous_read = read
                i = (i + 1) % 2

        except StopIteration:
            pass
