all: a.out merge.out
	rm -rf *.o *.mod
	mv a.out ../
	mv merge.out ../

#executable for dynamics
a.out: vbmode.o mapping.o main.o
	ifort -mkl -O2 -o a.out main.o mapping.o vbmode.o

#executable for merging
merge.out: merge.o
	ifort -o merge.out merge.o

#prerrequisites
mapping.o: mapping.f90
	ifort -c -O2 mapping.f90

vbmode.o: vbmode.f90
	ifort -c -O2 vbmode.f90

main.o: vbmode.o main.f90
	ifort -c -O2 main.f90

merge.o: merge.f90
	ifort -c -O2 merge.f90

#cleaning
clean:
	rm -rf *.o *.mod
