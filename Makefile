doc: qiimer-manual.pdf

qiimer-manual.pdf:
	R CMD roxygen -d .
	R CMD Rd2pdf --output=qiimer-manual.pdf --no-preview .

clean:
	rm -f qiimer-manual.pdf
	rm -f man/*.Rd
