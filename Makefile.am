# Copyright (C) 2011, 2012 Dominik Dahlem <dahlem@mit.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
ACLOCAL_AMFLAGS = -I m4 -I /usr/share/aclocal


SUBDIRS = m4 src

EXTRA_DIST = \
	README \
	$(srcdir)/m4 \
	$(srcdir)/mk


include mk/doxygen.am


DISTCHECK_CONFIGURE_FLAGS = --enable-valgrind

DISTCLEANFILES = configure
MAINTAINERCLEANFILES = "Makefile.in semantic.cache"
MOSTLYCLEANFILES = ${DX_CLEANFILES}


# add valgrind targets
#if VALGRIND
#valgrind: all
#	cd tests && $(MAKE) $(AM_MAKEFLAGS) valgrind
#valgrind.gen-suppressions: all
#	cd tests && $(MAKE) $(AM_MAKEFLAGS) valgrind.gen-suppressions
#else
#valgrind: ;
#valgrind.gen-suppressions: ;
#endif


.PHONY: valgrind valgrind.gen-suppressions \
	doxygen-run doxygen-doc $(DX_PS_GOAL) $(DX_PDF_GOAL)
