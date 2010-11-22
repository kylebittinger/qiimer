doc: ..pdf

..pdf:
	R CMD roxygen -d .
	R CMD Rd2pdf .

clean:
	rm -f ..pdf
	rm -f man/*.Rd
