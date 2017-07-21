# Makefile to create and install PaSDqc as a python package

.PHONY: dist

dist:
	python setup.py sdist
	pip install dist/*tar.gz
