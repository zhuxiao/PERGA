#!/bin/sh
	mkdir bin
	cd src/SRGA
	make
	mv -f SRGA ../../bin
	cd ../../src/SRGA_SF
	make
	mv -f SRGA_SF ../../bin
	cd ../../
	
