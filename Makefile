.PHONY: all for-commit reformat format isort rst unstaged todo mypy build pylint test tox docs
.PHONY: dist clean

.PHONY: tags flame sum

all: for-commit

for-commit: reformat format isort rst unstaged todo mypy pylint test docs tox

reformat: .clang-formatted
	autopep8 -i -r metacells tests setup.py

.clang-formatted: .clang-format metacells/extensions.cpp
	clang-format -i metacells/extensions.cpp
	@touch .clang-formatted

format:
	@if grep -n '^[ ].*[(%<>/+-]$$\|[^ =][[]$$' `git ls-files | grep '\.py$$'`; then false; fi
	@if grep -n '[^:`]`[^:`][^:`]*`[^`]' `git ls-files | grep '\.py$$' | grep -v metacells/scripts`; then false; fi

isort:
	isort metacells/**/*.py tests/*.py setup.py

unstaged:
	@if git status| grep 'not staged' > /dev/null; then git status; false; fi

todo:
	@if grep -i -n 'TODO''X\|reveal_''type' `git ls-files | grep -v pybind11`; then false; fi

mypy:
	mypy metacells tests

build:
	python setup.py build_ext --inplace

pylint: build
	pylint --jobs `nproc` metacells tests

test: build
	@rm -rf timing.*csv .coverage* coverage
	pytest -s -vv --cov=metacells tests

tox:
	tox

docs: sphinx/timing_script.rst
	sphinx-build -W -b html sphinx sphinx/generated

sphinx/timing_script.rst: \
    sphinx/timing_script.part.1 \
    sphinx/timing_script.part.2 \
    sphinx/timing_script.part.3 \
    sphinx/timing_script.part.4 \
    metacells/scripts/timing.py
	( cat sphinx/timing_script.part.1 \
	; python metacells/scripts/timing.py --help 2>&1 \
	| sed 's/timing.py/metacells_timing.py/;s/^/    /;s/`/``/g' \
	; cat sphinx/timing_script.part.2 \
	; python metacells/scripts/timing.py combine --help 2>&1 \
	| sed 's/timing.py/metacells_timing.py/;s/^/    /;s/`/``/g' \
	; cat sphinx/timing_script.part.3 \
	; python metacells/scripts/timing.py sum --help 2>&1 \
	| sed 's/timing.py/metacells_timing.py/;s/^/    /;s/`/``/g' \
	; cat sphinx/timing_script.part.4 \
	; python metacells/scripts/timing.py flame --help 2>&1 \
	| sed 's/timing.py/metacells_timing.py/;s/^/    /;s/`/``/g' \
	; cat sphinx/timing_script.part.5 \
	) > $@

rst: README.rst LICENSE.rst

LICENSE.rst: sphinx/license.rst
	cp $? $@

README.rst: prefix.rst sphinx/intro.rst sphinx/install.rst references.rst LICENSE.rst
	cat $? > $@

dist:
	python setup.py sdist

clean:
	rm -rf `cat .gitignore`

tags:
	ctags --python-kinds=-i -R metacells tests

coverage: coverage/index.html

coverage/index.html: timing.csv
	rm -rf coverage
	coverage html -d coverage

flame: flame.html

flame.html: timing.csv
	python metacells/scripts/timing.py combine timing.csv \
	| python metacells/scripts/timing.py flame -s \
	| flameview.py --sizename 'Total Elapsed Time' --sortby size \
	> flame.html

sum:
	python metacells/scripts/timing.py combine timing.csv \
	| python metacells/scripts/timing.py sum \
	| column -t -s,
