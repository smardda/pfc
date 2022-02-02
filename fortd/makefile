# config file.
include ../config/config.inc

src =  bsplvn.f bsplwn.f  bvalue.f  interu.f  interv.f  interw.f tb06a.f  tc01a.f  tg04a.f \
	isnumb.f spleval.f msrect.f lowor.f \
	isort.f dsort.f \
	cfftb1.f cffti1.f passb3.f passb5.f passf2.f passf4.f passf.f zfftf.f \
	cfftf1.f passb2.f passb4.f passb.f passf3.f passf5.f zfftb.f zffti.f


obj = $(src:.f=.o)

libsmarddalib$(DBG).a: $(obj) randlx.o
	ar ru libsmarddalib$(DBG).a $(obj) randlx.o
	ranlib libsmarddalib$(DBG).a

randlx.o:
	$(FC) $(FSFLAG) -c randlx.f

tar: $(src) README makefile
	tar cvf fort.tar $(src) randlx.f README makefile

clean: 
	rm -f $(obj) randlx.o

libclean: clean
	rm -f libsmarddalib$(DBG).a
