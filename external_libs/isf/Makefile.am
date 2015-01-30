if COMPILE_ISF
SUBDIRS = flib wrappers src
endif

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/report.py
