CRITERION_CFLAGS = $(shell pkg-config --cflags criterion)
CRITERION_LDLIBS = $(shell pkg-config --libs criterion)


ALL_TESTS := $(patsubst $(CURDIR)/%.c, %.out, $(wildcard $(CURDIR)/*.c))
CFLAGS = -g3 -O0 -Wno-cpp -lm

# se o criterion não largar a porta depois de ser finalizado, roda:
# fuser -k 1234/tcp
all: $(ALL_TESTS)
	$(foreach test, $^, ./$(test) $(ARGS);)

MEDIUM_TEST_DEP := $(CURDIR)/medium-tests.c $(PARENT_DIR)/src/medium.c $(PARENT_DIR)/src/includes/medium.h 
MEDIUM_TEST_DEP += $(PARENT_DIR)/src/argparse.c $(PARENT_DIR)/src/includes/argparse.h 
MEDIUM_TEST_DEP += $(PARENT_DIR)/src/includes/global.h $(PARENT_DIR)/src/includes/indexing.h $(PARENT_DIR)/src/includes/floatingpoint.h

medium-tests.out: $(MEDIUM_TEST_DEP)
	$(CC) -o $@ $^ $(CRITERION_CFLAGS) $(CRITERION_LDLIBS) $(STARPU_CFLAGS) $(STARPU_LDLIBS) -I $(PARENT_DIR)/$(INCLUDEDIR) $(CFLAGS)


KERNEL_TEST_DEP := $(CURDIR)/kernel-tests.c $(PARENT_DIR)/src/halo.c $(PARENT_DIR)/src/includes/halo.h 
KERNEL_TEST_DEP += $(PARENT_DIR)/src/kernel.c $(PARENT_DIR)/src/mem.c $(PARENT_DIR)/src/vector.c
KERNEL_TEST_DEP += $(PARENT_DIR)/src/medium.c $(PARENT_DIR)/src/argparse.c

kernel-tests.out: $(KERNEL_TEST_DEP)
	$(CC) -o $@ $^ $(CRITERION_CFLAGS) $(CRITERION_LDLIBS) $(STARPU_CFLAGS) $(STARPU_LDLIBS) -I $(PARENT_DIR)/$(INCLUDEDIR) $(CFLAGS)


HALO_TEST_DEP := $(CURDIR)/halo-tests.c $(PARENT_DIR)/src/halo.c $(PARENT_DIR)/src/includes/halo.h
HALO_TEST_DEP += $(PARENT_DIR)/src/includes/stencil.h $(PARENT_DIR)/src/includes/global.h $(PARENT_DIR)/src/includes/indexing.h

halo-tests.out: $(HALO_TEST_DEP)
	$(CC) -o $@ $^ $(CRITERION_CFLAGS) $(CRITERION_LDLIBS) -I $(PARENT_DIR)/$(INCLUDEDIR) $(CFLAGS)
