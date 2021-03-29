
.PHONY : batch combine download tables clean 

.SUFFIXES: .pdf

# --- GENERAL COMMANDS ---

clean :
	rm -rf output/*

# --- SERVER COMMANDS ---

# Submit sbatch script to slurm
batch :
	sbatch "code/batch/server.sbatch"

# Create raw tables from experiment results
# (run this after all experiments have run to completion)
combine :
	# mkdir -p output/tables
	sbatch "code/batch/combine_runs.sbatch"

# --- LOCAL COMMANDS ---

# Download results from server -- make sure to change spath as needed!
# (run this after running 'make combine' on the server)
spath := "$$MW:/home/livingstonb/GitHub/Discrete_HA/output/tables*"
cdate := $(shell date +"%m-%d-%Y-%T")
download :
	-mkdir -p output/server-$(cdate)
	-scp -r $(spath) output/server-$(cdate)/

# Construct latex tables and compile pdf from downloaded results
# (this code assumes that all of the downloaded table spreadsheets
# have been moved/copied to the 'output/tables' directory)
# tables :
# 	python code/+tables/create_tex_tables.py "output/tables"
# 	cp code/+tables/tables.tex output/tables/tables.tex
# 	cd output/tables && pdflatex tables
# 	cd output/tables && pdflatex tables

tables : tables1.pdf tables2.pdf tables3.pdf

%.pdf :
	python code/+tables/create_tex_tables.py "output/tables/$*"
	cp code/+tables/tables.tex output/tables/$*/$*.tex
	cd output/tables/$* && pdflatex $*
	cd output/tables/$* && pdflatex $*