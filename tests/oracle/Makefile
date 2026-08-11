# Differential oracle for the expression parser.
#
#   make CFITSIO=/path/to/cfitsio        build ./oracle
#   make CFITSIO=/path/to/cfitsio run     dump the corpus through the C parser
#   make CFITSIO=/path/to/cfitsio check   diff the C parser against our golden
#   make CFITSIO=/path/to/cfitsio lexdump build the token-level oracle
#
# CFITSIO must point at a *built* CFITSIO source tree (one where `make` has
# produced .libs/libcfitsio.so).  See README.md.

CFITSIO ?= $(HOME)/code/cfitsio
LIBDIR   = $(CFITSIO)/.libs
LIB      = $(firstword $(wildcard $(LIBDIR)/libcfitsio.so.* $(LIBDIR)/libcfitsio.so))

CORPUS   = ../fixtures/eval_corpus.txt
GOLDEN   = ../fixtures/eval_corpus.golden

CFLAGS  += -O2 -Wall -I$(CFITSIO)
LDFLAGS += -lm -Wl,-rpath,$(LIBDIR)

.PHONY: all run check lexdump clean

all: oracle

oracle: oracle.c
	@test -n "$(LIB)" || { echo "no libcfitsio in $(LIBDIR); build CFITSIO first"; exit 1; }
	$(CC) $(CFLAGS) -o $@ $< $(LIB) $(LDFLAGS)

run: oracle
	./oracle $(CORPUS)

# Lines where the C library and rsfitsio are known to disagree; see README.md.
check: oracle
	@./oracle $(CORPUS) > /tmp/oracle.out
	@diff -u $(GOLDEN) /tmp/oracle.out > /tmp/oracle.diff \
		&& echo "identical to $(GOLDEN)" \
		|| { echo "differences (see README.md for the known ones):"; \
		     cat /tmp/oracle.diff; }

lexdump: lexdump.l
	flex -o lexdump.c lexdump.l
	$(CC) -O2 -o lexdump lexdump.c

clean:
	rm -f oracle lexdump lexdump.c
