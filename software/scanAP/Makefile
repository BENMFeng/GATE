PRGS = scanAP
CC = g++
CFLAGS = -O3

build: $(PRGS)

scanAP: scanAP.cpp local_alignment.cpp
	$(CC) $(CFLAGS) -o $@ $^

cp:
	cp $(PRGS) ../../bin/
	cp ./*.fa ../../bin/
	cp ./align.mat ../../bin/

install: build cp

clean:
	rm $(PRGS)
