CURSO = 25-26
Calendario-Econometria-Aplicada.pdf: $(CURSO)/Calendario-Econometria-Aplicada.tex $(CURSO)/CalendarioEscolar.tex $(CURSO)/clases.tex
	cd $(CURSO) && pdflatex -interaction batchmode Calendario-Econometria-Aplicada.tex
	cd $(CURSO) && rm -f Calendario-Econometria-Aplicada.log
	cd $(CURSO) && rm -f Calendario-Econometria-Aplicada.aux
	mv $(CURSO)/Calendario-Econometria-Aplicada.pdf ../docs
