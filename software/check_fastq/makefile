PRGS = checkfastq
CC = g++
CFLAGS = -O3

build: $(PRGS)

checkfastq: checkfastq.cpp
	$(CC) $(CFLAGS) -o $@   $^

cp:
	cp $(PRGS) ../bin/
	cp distribute_fqcheck.pl ../bin/

install: build cp

clean:
	rm $(PRGS)
