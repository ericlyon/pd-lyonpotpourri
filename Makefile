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
	chameleon~.c \
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

## sources for the shared library
# fft.c       :?: convert.c convolver~.c leanunconvert.c     splitbank~.c
# fft4.c      :?: buffet~.c convolver~.c magfreq_analysis~.c splitbank~.c
# fold.c        : buffet~.c              magfreq_analysis~.c splitbank~.c !squash~.c
# convert.c     : buffet~.c              magfreq_analysis~.c !splitbank~.c
# makewindows.c : buffet~.c magfreq_analysis~.c splitbank~.c !squash~.c
# power_of_two.c:           magfreq_analysis~.c
# from_msp.c: adsr~.c bashfest~.c buffet~.c counter~.c stutter~.c vdp~.c

shared.sources = \
	convert.c \
	fft.c \
	fft4.c \
	fold.c \
	makewindows.c \
	power_of_two.c \
	from_msp.c \
	$(empty)

unused_shared_sources = \
	PenroseRand.c \
	PenroseOscil.c \
	bloscbank.c \
	overlapadd.c \
	leanconvert.c \
	leanunconvert.c \
	unconvert.c \
	qsortE.c \
	$(empty)

# extra files
datafiles = \
	$(wildcard *.pd) \
	CHANGELOG.txt \
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
