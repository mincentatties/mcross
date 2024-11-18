# Makefile for McRoss for mac
#----------------------------

FFPOBJECTS = ffp.o formfactor.o ff_table.o
NUCOBJECTS = mcross.o nucsro.o
MAGOBJECTS = mcross.o magsro.o formfactor.o ff_table.o
WCOBJ      = mcross.o wc_calc.o
COROBJ     = mcross.o cor_calc.o

# debugging compiler command
# COML = gfortran -w -v -c -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic 
# optimized code command
COML = gfortran -w -v -c -O3#-ftree-vectorize 
LINK = gfortran -o

all:	ffp nucsro magsro wc_calc cor_calc

ffp:	$(FFPOBJECTS)
		$(LINK) ffp/ffp $(FFPOBJECTS)
nucsro:	$(NUCOBJECTS)
		$(LINK) nucsro/nucsro $(NUCOBJECTS)
magsro:	$(MAGOBJECTS)
		$(LINK) magsro/magsro $(MAGOBJECTS)	
wc_calc:$(WCOBJ)
	    $(LINK) nucsro/wc_calc $(WCOBJ)
cor_calc:$(COROBJ)
	    $(LINK) magsro/cor_calc $(COROBJ)


ffp.o: 			ffp/ffp.f90
			$(COML) ffp/ffp.f90

formfactor.o: 		ffp/formfactor.f90
			$(COML) ffp/formfactor.f90

ff_table.o: 		ffp/ff_table.f90
			$(COML) ffp/ff_table.f90

nucsro.o:		nucsro/nucsro.f90
			$(COML) nucsro/nucsro.f90
		
magsro.o:		magsro/magsro.f90
			$(COML) magsro/magsro.f90
		
mcross.o:		mcross.f90
			$(COML) mcross.f90

wc_calc.o:		nucsro/wc_calc.f90
			$(COML) nucsro/wc_calc.f90

cor_calc.o:		magsro/cor_calc.f90
			$(COML) magsro/cor_calc.f90

clean: 
		rm -f $(FFPOBJECTS)
		rm -f $(NUCOBJECTS)
		rm -f $(MAGOBJECTS)
		rm -f $(WCOBJ)
		rm -f $(COROBJ)
		rm -f mcross.mod ffp/ffp magsro/magsro nucsro/nucsro magsro/cor_calc nucsro/wc_calc
		rm -f ffp/ffp magsro/magsro nucsro/nucsro nucsro/wc_calc magsro/cor_calc

install:
	@echo Must be administrator to install programs
		sudo rm /usr/local/bin/ffp
		sudo rm /usr/local/bin/magsro
		sudo rm /usr/local/bin/nucsro
		sudo rm /usr/local/bin/wc_calc
		sudo rm /usr/local/bin/cor_calc
		sudo cp ffp/ffp       	/usr/local/bin/ffp
		sudo cp magsro/magsro 	/usr/local/bin/magsro
		sudo cp nucsro/nucsro 	/usr/local/bin/nucsro
		sudo cp nucsro/wc_calc  /usr/local/bin/wc_calc
		sudo cp magsro/cor_calc /usr/local/bin/cor_calc

		
