#/usr/bin/make

GLOTZER_HOME:=$$WORK/glotzer

qm: QuaternionMultiply.cc
	c++ QuaternionMultiply.cc -I$(GLOTZER_HOME)/hoomd-plugin-hpmc/cppmodule -I$(GLOTZER_HOME)/hoomd-install/include/ -I$(TACC_CUDA_INC) -I/opt/apps/intel13/mvapich2/1.9/include \
		-I$(BOOST_ROOT)/include -o qm

pos2pov: pos2pov.cc
	c++ pos2pov.cc -I$(GLOTZER_HOME)/hoomd-plugin-hpmc/cppmodule -I$(GLOTZER_HOME)/hoomd-install/include/ -I$(TACC_CUDA_INC) -I/opt/apps/intel13/mvapich2/1.9/include \
		-I$(BOOST_ROOT)/include -o pos2pov

echo:	
	@echo GLOTZER_HOME=$(GLOTZER_HOME)
	@echo TACC_CUDA_INC=$(TACC_CUDA_INC)
	@echo BOOST_ROOT=$(BOOST_ROOT)
