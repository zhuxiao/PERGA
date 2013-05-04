#!/bin/sh
	mkdir bin
	cd src/SRGA
	make
	mv -f srga ../../bin
	cp -r model ../../bin
	cd ../../src/SRGA_SF
	make
	mv -f srga_sf ../../bin
	cd ../../
	
