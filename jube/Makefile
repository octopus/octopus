all: bench

clean:
	rm -rf *~ xmllogs logs benchlog tmp results .bench_current_id.dat svn-commit.*

clean_run: clean
	rm -rf run

code:
	@if [ -d octopus ]; then \
	  if [ -f octopus/Makefile ]; then \
	    (cd octopus; make clean); \
	  fi; \
	  (cd octopus; svn update); \
	else \
	  svn co http://www.tddft.org/svn/octopus/trunk octopus; \
	fi

bench: code
	case `hostname` in \
	  jugene* ) platform='IBM-BGP-Jugene'   ;; \
	  louhi*  ) platform='Cray-XT4-Louhi'   ;; \
	  p6012*  ) platform='IBM-SP6-Huygens'  ;; \
	  s1c1b*  ) platform='IBM-Cell-MariCel' ;; \
	esac; \
	if [ -z "$$platform" ]; then \
	  echo "Platform '`hostname`' not found"; exit 1; \
        else \
	  echo Platform is: $$platform; \
	  ../../bench/jube ./nbench-$$platform.xml; \
	fi;

