#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Writing trees in Newick format.

"""

__docformat__ = 'restructuredtext'


### IMPORTS

import re

from phylo.core import io

__all__ = [
	'WriterDialect',
]


### CONSTANTS & DEFINES

### IMPLEMENTATION ####


class WriterDialect (io.Dialect):
	defaults = {
		'support_format': '%.2f',
		'distance_format': '%.3f',
		'allow_singletons': False,
		'include_distances': True,
		'include_support': True,
		'quote_all_names': False,
	}
		
	def validate (self):
		pass

