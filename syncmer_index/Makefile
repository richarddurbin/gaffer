all: syncmer_index

gsacak.o: gsa-is/gsacak.c
	g++ -O3 -c $<

malloc_count.o: gsa-is/experiments/external/malloc_count/malloc_count.c
	g++ -O3 -c $<

syncmer_index: syncmer_index.cpp gsacak.o malloc_count.o syncmer_index.cpp
	g++ -O3 -I./gsa-is -o syncmer_index gsacak.o malloc_count.o ../ONElib.o ../utils.o syncmer_index.cpp -ldl -lz

clean:
	rm gsacak.o malloc_count.o syncmer_index

#421  g++ -c gsacak.c experiments/external/malloc_count/malloc_count.c
#421  g++ -c gsacak.c experiments/external/malloc_count/malloc_count.c
#  422  g++ test2.cpp -o test2 gsacak.o malloc_count.o -ldl
