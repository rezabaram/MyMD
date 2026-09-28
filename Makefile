include Makefile.inc
DIRS= tools


run: ellipmd
	time bin/run.sh $(CONFIG)

ellipmd:	*.cc include/*.h 
	$(CC)  main.cc  $(FLAGS) $(WARNFLAGS) $(OPTFLAGS) -o ellipmd $(LDFLAGS)

# Physics regression check -- run this before and after any change to the
# solver.  See bench/check_physics.py and ROADMAP.md (Phase 0).
check: ellipmd
	python3 bench/check_physics.py

# Sanitizer build.  Uses clang rather than the default compiler because the
# Homebrew GCC on this machine ships no sanitizer runtime, and clang could not
# build this code at all before the C++17 migration.
#
# This build found two real bugs on its first run: a heap-buffer-overflow in
# DisBetaDistribution (bins was allocated with nbins slots but every loop
# indexed 0..nbins inclusive) and an invalid downcast in
# CBaseConfig::get_param<T>.  Run it after touching anything that manages
# memory or does a cast.
SAN_CC     ?= clang++
SAN_FLAGS  ?= -std=c++17 -g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer
asan: main.cc
	$(SAN_CC) main.cc $(FLAGS) -o ellipmd-asan $(LDFLAGS) $(SAN_FLAGS)
	@echo "built ./ellipmd-asan -- run:  ./ellipmd-asan <seed> <config>"

# Same, but checking the reference cases automatically.
asan-check: asan
	@for c in $$(ls bench/reference); do \
		d=$$(mktemp -d); \
		cp bench/reference/$$c/* $$d/ 2>/dev/null; \
		( cd $$d && UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0 \
		  ASAN_OPTIONS=abort_on_error=0 \
		  $(ROOT)/ellipmd-asan 1 config > out.log 2>&1 ); \
		n=$$(grep -cE 'ERROR: AddressSanitizer|runtime error:' $$d/out.log); \
		printf '  %-16s findings=%s\n' $$c $$n; \
		rm -rf $$d; \
	done

# Timing benchmark.  Sequential on purpose; writes bench/results.json.
#   make bench              # all four configurations (~16 min)
#   make bench BENCH_ARGS="--only B2"    # ~80 s
BENCH_ARGS ?=
bench: ellipmd
	python3 bench/run_bench.py $(BENCH_ARGS)


update:
	git pull origin master
.PHONY: all run clean update tools gsl deps viewer viewer-check ovito ovito-render dump check bench asan asan-check live viewer-smoke
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

# Live dashboard: edit the parameters, press Start, and watch the run fill in
# the browser.  Local only -- the server is the Python standard library.
#   make live                        # opens on http://127.0.0.1:8770/
#   make live LIVE_CONFIG=config_quick PORT=8800
LIVE_CONFIG ?= $(ROOT)/config_rain
PORT ?= 8770
live:
	python3 tools/live_viewer.py --port $(PORT) --config $(LIVE_CONFIG)

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

# Execute the generated viewer pages under node with a stubbed browser.  This
# is the only check that runs the viewer JavaScript rather than just parsing
# it, and it exists because two runtime errors reached a browser from here
# (a temporal-dead-zone ReferenceError and an undeclared variable) that
# `node --check` cannot see.  Skipped if node is not installed.
viewer-smoke:
	python3 tools/viewer_smoketest.py

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
