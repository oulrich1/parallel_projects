# GNU Make solution makefile autogenerated by Premake
# Type "make help" for usage help

ifndef config
  config=release
endif
export config

PROJECTS := matrixmult

.PHONY: all clean help $(PROJECTS)

all: $(PROJECTS)

matrixmult: 
	@echo "==== Building matrixmult ($(config)) ===="
	@${MAKE} --no-print-directory -C build -f Makefile

clean:
	@${MAKE} --no-print-directory -C build -f Makefile clean

help:
	@echo "Usage: make [config=name] [target]"
	@echo ""
	@echo "CONFIGURATIONS:"
	@echo "   release"
	@echo "   debug"
	@echo ""
	@echo "TARGETS:"
	@echo "   all (default)"
	@echo "   clean"
	@echo "   matrixmult"
	@echo ""
	@echo "For more information, see http://industriousone.com/premake/quick-start"
