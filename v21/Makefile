all: build/v21.pdf

# hier Python-Skripte:
build/Isotop1.pdf build/Isotop2.pdf: plot.py ../matplotlibrc ../header-matplotlib.tex | build
	# so that matplotlib can find the tex header when running
	# LaTeX in the tmp directory
	# and set the matplotlibrc
	TEXINPUTS=$$(pwd)/..: MATPLOTLIBRC=../matplotlibrc python plot.py

# hier weitere Abhängigkeiten für build/v21.pdf deklarieren:
build/v21.pdf: build/Isotop1.pdf build/Isotop2.pdf

build/v21.pdf: FORCE | build
	# to find header and bib files in the main directory
	TEXINPUTS=..: \
	BIBINPUTS=..: \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	v21.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
