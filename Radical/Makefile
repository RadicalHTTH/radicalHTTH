# Compiler and flags
CC = gcc
CFLAGS = -Wall -g -Wno-empty-body -Wno-unused-but-set-variable -Wno-unused-variable -Wno-address

# Directory structure
SRCDIR = src
OBJDIR = obj
BINDIR = bin

# Source and object files
SRC = $(SRCDIR)/evoradical_re.c $(SRCDIR)/eigen.c $(SRCDIR)/toolsfromPAML.c
OBJ = $(OBJDIR)/evoradical_re.o $(OBJDIR)/eigen.o $(OBJDIR)/toolsfromPAML.o

# Output File
TARGET = $(BINDIR)/reEvoRadical

# make rules
all: $(BINDIR) $(OBJDIR) $(TARGET)

# Generating an executable file
$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ -lm

# Generate .o files (include xdfpmin.c in dependencies)
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(SRCDIR)/*.h $(SRCDIR)/xdfpmin.c | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Create the directory if it does not exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

# Clean up
clean:
	rm -rf $(OBJDIR) $(BINDIR)
