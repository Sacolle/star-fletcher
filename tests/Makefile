CRITERION_CFLAGS = $(shell pkg-config --cflags criterion)
CRITERION_LDLIBS = $(shell pkg-config --libs criterion)


ALL_TESTS := $(patsubst $(CURDIR)/%.c, %.out, $(wildcard $(CURDIR)/*.c))
CFLAGS = -g3 -O0 -Wno-cpp -lm

# Writes <test>.d with the project headers included by every source of the test, so that editing
# a header rebuilds the tests that use it. Each test is compiled and linked in a single command,
# where -MMD can't be used (each source would overwrite the previous .d), hence the -MM pass.
# Sources compiled into a test: the .c prerequisites, minus the generated cross-deriv.gen.c
# (listed by the .d files, but only meant to be included by derivatives-impl.h).
TEST_SRCS = $(filter-out %.gen.c,$(filter %.c,$^))

DEPGEN = $(CC) -MM -MP -MT $@ $(TEST_SRCS) -I $(PARENT_DIR)/$(INCLUDEDIR) $(CRITERION_CFLAGS) $(STARPU_CFLAGS) > $(@:.out=.d)

.PHONY: all clean

# se o criterion não largar a porta depois de ser finalizado, roda:
# fuser -k 1234/tcp
all: $(ALL_TESTS)
	$(foreach test, $^, ./$(test) $(ARGS);)

deriv-tests.out: $(CURDIR)/deriv-tests.c $(PARENT_DIR)/src/derivatives.c $(PARENT_DIR)/src/includes/derivatives.h 
	@$(DEPGEN)
	$(CC) -o $@ $(TEST_SRCS) $(CRITERION_CFLAGS) $(CRITERION_LDLIBS) -I $(PARENT_DIR)/$(INCLUDEDIR) $(CFLAGS)


MEDIUM_TEST_DEP := $(CURDIR)/medium-tests.c $(PARENT_DIR)/src/medium.c $(PARENT_DIR)/src/includes/medium.h 
MEDIUM_TEST_DEP += $(PARENT_DIR)/src/argparse.c $(PARENT_DIR)/src/includes/argparse.h 
MEDIUM_TEST_DEP += $(PARENT_DIR)/src/includes/macros.h $(PARENT_DIR)/src/includes/floatingpoint.h

medium-tests.out: $(MEDIUM_TEST_DEP)
	@$(DEPGEN)
	$(CC) -o $@ $(TEST_SRCS) $(CRITERION_CFLAGS) $(CRITERION_LDLIBS) $(STARPU_CFLAGS) $(STARPU_LDLIBS) -I $(PARENT_DIR)/$(INCLUDEDIR) $(CFLAGS)


KERNEL_TEST_DEP := $(CURDIR)/kernel-tests.c $(PARENT_DIR)/src/derivatives.c $(PARENT_DIR)/src/includes/derivatives.h 
KERNEL_TEST_DEP += $(PARENT_DIR)/src/kernel.c $(PARENT_DIR)/src/mem.c $(PARENT_DIR)/src/vector.c
KERNEL_TEST_DEP += $(PARENT_DIR)/src/medium.c $(PARENT_DIR)/src/argparse.c

kernel-tests.out: $(KERNEL_TEST_DEP)
	@$(DEPGEN)
	$(CC) -o $@ $(TEST_SRCS) $(CRITERION_CFLAGS) $(CRITERION_LDLIBS) $(STARPU_CFLAGS) $(STARPU_LDLIBS) -I $(PARENT_DIR)/$(INCLUDEDIR) $(CFLAGS)

clean:
	rm -f $(ALL_TESTS) $(ALL_TESTS:.out=.d)

# header dependencies written by DEPGEN (missing on the first build, which compiles everything anyway)
-include $(ALL_TESTS:.out=.d)
