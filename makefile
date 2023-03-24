REP=$(shell basename $(PWD))
SOURCES=main.c
CIBLE= TP5.exe
CFLAGS=-Wall -Wno-format-overflow

# makefile générique pour produire un code source 
# dont le nom correspond au nom du répertoire qui le contient

all:
	@echo -n "Production de $(CIBLE)"
	@echo " à partir des fichiers : $(SOURCES)"
	gcc $(CFLAGS) $(SOURCES) -o $(CIBLE)
	@echo "Le programme $(CIBLE) a été produit dans le répertoire $(REP)"


clean: 
	@echo "Nettoyage de $(CIBLE)"
	@rm -rf $(CIBLE)
