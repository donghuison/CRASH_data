# Makefile for CRASH opacity test
# Compiler settings
FC = gfortran
FFLAGS = -O2 -Wall -fcheck=all -g
LDFLAGS = 

# Target executable
TARGET = test_crash_opacity

# Source files
SOURCES = crash_opacity_reader.f90 test_crash_opacity.f90
OBJECTS = crash_opacity_reader.o test_crash_opacity.o

# Default target
all: $(TARGET)

# Build the executable
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)
	@echo "Build complete: $(TARGET)"

# Compile module
crash_opacity_reader.o: crash_opacity_reader.f90
	$(FC) $(FFLAGS) -c $<

# Compile main program (depends on module)
test_crash_opacity.o: test_crash_opacity.f90 crash_opacity_reader.o
	$(FC) $(FFLAGS) -c $<

# Run the test
run: $(TARGET)
	@echo "Running CRASH opacity test..."
	./$(TARGET)

# Clean build files
clean:
	rm -f *.o *.mod $(TARGET)
	@echo "Clean complete"

# Rebuild everything
rebuild: clean all

.PHONY: all clean run rebuild