#/usr/bin/make

GLOTZER_HOME:=$$WORK/glotzer
VACUUMMS_HOME:=$$HOME

qm: QuaternionMultiply.cc
	c++ QuaternionMultiply.cc -I$(GLOTZER_HOME)/hoomd-plugin-hpmc/cppmodule -I$(GLOTZER_HOME)/hoomd-install/include/ -I$(TACC_CUDA_INC) -I/opt/apps/intel13/mvapich2/1.9/include \
		-I$(BOOST_ROOT)/include -o qm

pos2pov: pos2pov.cc
	c++ pos2pov.cc -I$(GLOTZER_HOME)/hoomd-plugin-hpmc/cppmodule -I$(GLOTZER_HOME)/hoomd-install/include/ -I$(TACC_CUDA_INC) -I/opt/apps/intel13/mvapich2/1.9/include \
		-I$(BOOST_ROOT)/include -o pos2pov

pmft2pov: pmft2pov.c
	gcc -g pmft2pov.c \
		-I$(GLOTZER_HOME)/hoomd-plugin-hpmc/cppmodule -I$(GLOTZER_HOME)/hoomd-install/include/ \
                -I$(TACC_CUDA_INC) -I/opt/apps/intel13/mvapich2/1.9/include \
                -I$(VACUUMMS_HOME)/usr/include -L$(VACUUMMS_HOME)/usr/lib \
		-lm -lftw -lftw_pov\
		-o pmft2pov

rattle: rattle.cc
	c++ -g rattle.cc -I$(GLOTZER_HOME)/hoomd-plugin-hpmc/cppmodule -I$(GLOTZER_HOME)/hoomd-install/include/ \
                -I$(TACC_CUDA_INC) -I/opt/apps/intel13/mvapich2/1.9/include \
                -I$(VACUUMMS_HOME)/usr/include -L$(VACUUMMS_HOME)/usr/lib \
		-I$(BOOST_ROOT)/include -o rattle

testit: testit.cc
	c++ testit.cc -I$(GLOTZER_HOME)/hoomd-plugin-hpmc/cppmodule -I$(GLOTZER_HOME)/hoomd-install/include/ -I$(TACC_CUDA_INC) -I/opt/apps/intel13/mvapich2/1.9/include \
		-I$(BOOST_ROOT)/include -o testit

#install: pos2pov pmft2pov rattle
install: pos2pov pmft2pov 
	cp pmft2pov pos2pov ../bin/
#	cp pmft2pov pos2pov rattle ../bin/

echo:	
	@echo GLOTZER_HOME=$(GLOTZER_HOME)
	@echo TACC_CUDA_INC=$(TACC_CUDA_INC)
	@echo BOOST_ROOT=$(BOOST_ROOT)
