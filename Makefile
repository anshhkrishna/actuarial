.PHONY: all install download verify build analysis clean

PY := python

all:
	$(PY) run_all.py

install:
	$(PY) -m pip install -r requirements.txt

download:
	$(PY) src/download.py

verify:
	$(PY) src/verify_variables.py

build:
	$(PY) src/build_cohort.py

analysis:
	$(PY) src/analysis.py

clean:
	rm -f data/processed/*.csv data/processed/*.json
	rm -f results/*.csv results/*.txt
	rm -f figures/*.png
