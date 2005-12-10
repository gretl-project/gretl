# location of various things

WEBDIR = /home/cottrell/stats/esl/website

REF_CHAPTERS = chapters/cmdref.xml chapters/optarg.xml
CMDLIST = chapters/cmdlist.xml
CMDREF = commands/gretl_commands.xml

REF_CHAPTERS_IT = chapters_it/cmdref.xml chapters_it/optarg.xml
CMDLIST_IT = chapters_it/cmdlist.xml
CMDREF_IT = commands/gretl_commands.xml

SUBDIRS = tex tex_it

.PHONY : subdirs $(SUBDIRS) texfigs clean install 

subdirs: $(SUBDIRS)

$(SUBDIRS):
	make -C $@

commands/eqn_list.xml:
	make -C commands -f cmdref.mk eqn_list.xml	

texfigs: commands/eqn_list.xml
	./texmath2png.pl commands/eqn_list.xml

clean:
	make -C tex clean
	make -C tex_it clean
	make -C tex_es clean

install: 
	make -C tex install
	make -C tex_it install
	make -C tex_es clean
