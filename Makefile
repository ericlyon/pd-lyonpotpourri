# Makefile to build the 'lyonpotpourri' library for Pure Data.
# Needs Makefile.pdlibbuilder as helper makefile for platform-dependent build
# settings and rules.

# library name
lib.name = lyonpotpourri

# sources for the objectclasses
class.sources = \
	adsr~.c \
	buffet~.c \
	bvplay~.c \
	channel~.c \
	chopper~.c \
	clean_selector~.c \
	click~.c \
	click2bang~.c \
	click2float~.c \
	clickhold~.c \
	distortion~.c \
	dmach~.c \
	expflam~.c \
	flanjah~.c \
	function~.c \
	granola~.c \
	granulesf~.c \
	granule~.c \
	kbuffer~.c \
	killdc~.c \
	magfreq_analysis~.c \
	markov~.c \
	mask~.c \
	oscil~.c \
	phasemod~.c \
	player~.c \
	pulser~.c \
	rtrig~.c \
	samm~.c \
	sigseq~.c \
	vdb~.c \
	vdp~.c \
	waveshape~.c \
	epluribus~.c \
	dynss~.c \
	counter~.c \
	latch~.c \
	sarec~.c \
	convolver~.c \
	npan~.c \
	shoehorn~.c \
	rotapan~.c \
	sel~.c \
	squash~.c \
	windowvec~.c \
	cartopol~.c \
	poltocar~.c \
	arrayfilt~.c \
	splitspec~.c \
	stutter~.c \
	vecdex~.c \
	quadpan~.c \
	splitbank~.c \
	$(empty)

bashfest~.class.sources = \
	bashfest~.c \
	bashfest_dsp.c \
	bashfest_helper.c \
	ellipse.c \
	$(empty)

# sources for the shared library
shared.sources = \
	PenroseOscil.c \
	PenroseRand.c \
	bloscbank.c \
	convert.c \
	fft.c \
	fft4.c \
	fold.c \
	leanconvert.c \
	leanunconvert.c \
	makewindows.c \
	overlapadd.c \
	power_of_two.c \
	qsortE.c \
	unconvert.c \
	from_msp.c \
	$(empty)

# extra files
datafiles = \
	$(wildcard *.pd) \
	LICENSE.txt \
	README.txt \
	lyonpotpourri-meta.pd \
	$(empty)
# extra dirs
datadirs = examples \
	$(empty)

##############################################################################################
# include the actual build-system
PDLIBBUILDER_DIR=pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
