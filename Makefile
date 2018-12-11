CC=gcc -O3 -Wall

main: src/*.c
	$(CC) -lm src/*.c -o main

.PHONY: clean
clean:
	rm main

