batch :
	sbatch "code/batch/server.sbatch"

combine :
	mkdir -p output/tables
	sbatch "code/batch/combine_runs.sbatch"

spath := "$$MW:/home/livingstonb/GitHub/Discrete_HA/output/tables/*"
cdate := $(shell date +"%m-%d-%Y-%T")
download :
	-mkdir -p output/server-$(cdate)
	-scp $(spath) output/server-$(cdate)

tables :
	python code/+tables/create_tex_tables.py "output/tables"
	cp code/+tables/tables.tex output/tables/tables.tex
	cd output/tables && pdflatex tables
	cd output/tables && pdflatex tables

clean :
	rm -rf output/*

.PHONY : batch combine download tables clean 