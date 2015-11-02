#TYPE =	float
TYPE =	double

testfpzip: testfpzip.cpp ../lib/libfpzip.a
	$(CXX) -ansi -Wall -g -DTYPE=$(TYPE) testfpzip.cpp -I../inc -L../lib -lfpzip -o testfpzip

clean:
	@rm -f testfpzip
