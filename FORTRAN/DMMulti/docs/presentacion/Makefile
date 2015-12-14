
slides.pdf: gnuplot src/slides.tex
	mkdir -p build
	cd build; pdflatex --shell-escape ../src/slides.tex; pdflatex --shell-escape ../src/slides.tex; cp slides.pdf ..

clean:
	rm -rf build; rm -f slides.pdf

distclean: clean
	cd src/src_figs; 

neat:
	rm -rf build

.PHONY: gnuplot clean distclean neat
