SRGA: srgaMain.o srga.o graph.o contig.o update.o  util.o hashPE.o PEAssembly.o
	gcc -g -o $@ $^ -lm  # $@:目标文件; $^:所有依赖文件
	rm *.o
.c.o:  #所有的.o文件都要依赖于相应的.c文件
	gcc -g -c $< -lm

clean:
	rm *.o
	rm *~
	rm *.bak
