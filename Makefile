include Makefile.inc
DIRS= tools


run: ellipmd
	time bin/run.sh $(CONFIG)

ellipmd:	*.cc include/*.h 
	$(CC)  main.cc  $(FLAGS) $(DEBUGFLAGS) -o ellipmd $(LDFLAGS)

# Physics regression check -- run this before and after any change to the
# solver.  See bench/check_physics.py and ROADMAP.md (Phase 0).
check: ellipmd
	python3 bench/check_physics.py

# Timing benchmark.  Sequential on purpose; writes bench/results.json.
#   make bench              # all four configurations (~16 min)
#   make bench BENCH_ARGS="--only B2"    # ~80 s
BENCH_ARGS ?=
bench: ellipmd
	python3 bench/run_bench.py $(BENCH_ARGS)


update:
	git pull origin master
.PHONY: all run clean update tools gsl deps viewer viewer-check ovito ovito-render dump check bench
tools:
	$(MAKE) -C tools

# ---- modern visualisation (see VISUALIZATION.md) -------------------------
# Self-contained interactive HTML viewer, no installs needed:
#   make viewer OUT='out0*'
OUT ?= out0* outend
HTML ?= trajectory.html
COLOR ?= uniform
viewer:
	python3 tools/web_viewer.py $(OUT) -o $(HTML) --color $(COLOR)
	@echo "open file://$$(cd $$(dirname '$(HTML)') && pwd)/$$(basename '$(HTML)')"

# LAMMPS-dump conversion, for dragging into the OVITO GUI:
#   make dump OUT='out0*' DUMP=trajectory.dump
DUMP ?= trajectory.dump
dump:
	python3 tools/snapshot_to_dump.py $(OUT) -o $(DUMP)
	@echo "open it with:  open -a Ovito $(DUMP)    (brew install --cask ovito)"

# Known-answer page: ellipsoids whose long axes should point along +X, +Y, -Z,
# then a triaxial one and a c/a ladder.  Use it to confirm the viewer is not
# lying about orientation.
viewer-check:
	python3 tools/orientation_sample.py "$(ROOT)/data/orientation_sample"
	python3 tools/web_viewer.py "$(ROOT)/data/orientation_sample" \
		-o orientation_sample.html \
		--title "Orientation check (long axes: X, Y, Z, triaxial, then c/a = 1,2,3)"
	@echo "open file://$(ROOT)/orientation_sample.html"

# OVITO (https://ovito.org) python module in a local virtualenv:
#   make ovito
#   make ovito-render FILE=out00010 PNG=frame.png
OVITO ?= $(ROOT)/.deps/venv/bin/python
ovito:
	@mkdir -p "$(ROOT)/.deps"
	@test -x "$(ROOT)/.deps/venv/bin/python" || \
		UV_CACHE_DIR="$(ROOT)/.deps/uvcache" uv venv "$(ROOT)/.deps/venv"
	UV_CACHE_DIR="$(ROOT)/.deps/uvcache" uv pip install --python "$(OVITO)" ovito

FILE ?= out00010
PNG ?= frame.png
OVITO_FLAGS ?=
ovito-render: ovito
	$(OVITO) tools/ovito_reader.py $(FILE) --out $(PNG) $(OVITO_FLAGS)

# Fetch and build the GSL dependency into .deps/gsl (used automatically by
# Makefile.inc).  Skip this if GSL is already installed system-wide, e.g.
# 'brew install gsl', in which case the Homebrew copy is picked up instead.
gsl deps:
	@if [ -d "$(GSL_PREFIX)/include/gsl" ]; then \
		echo "GSL already built in $(GSL_PREFIX)"; \
	else \
		set -e; \
		mkdir -p "$(ROOT)/.deps"; \
		cd "$(ROOT)/.deps"; \
		echo "Downloading GSL 2.8 ..."; \
		curl -sSL -o gsl.tar.gz https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz; \
		tar xzf gsl.tar.gz; \
		cd gsl-2.8; \
		./configure --prefix="$(GSL_PREFIX)" --disable-shared --enable-static > configure.log 2>&1; \
		make -j8 > build.log 2>&1; \
		make install > install.log 2>&1; \
		echo "GSL installed in $(GSL_PREFIX)"; \
	fi

clean:
	rm -f log_energy log out0* outend *jpg test.avi

zip:
	zip md.zip *.cc include/*.h Makefile Makefile.inc bin/run.sh config

force_look :
	true
