TESTS=test_0 test_1 test_2
all: ${TESTS}

test_%: reproducer.cc
	mpicxx $^ -o $@ -DTEST_CASE=${patsubst test_%,%,$@} \
		-Wl,-rpath,${shell pkg-config --variable=libdir yaxt} \
		${shell pkg-config --cflags --libs yac} \
		${shell pkg-config --cflags --libs yaxt} \
		${shell pkg-config --cflags --libs proj}

clean:
	@rm -f ${TESTS}
